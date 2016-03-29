package com.lifetechnologies.solid.wt.mapper;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.IndexedFastaDatabase;
import com.lifetechnologies.solid.wt.ReadMapper;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.splice.SpliceJunctionExtractor;
import com.lifetechnologies.util.FileHeader;
import com.lifetechnologies.util.FileUtils;
import com.lifetechnologies.util.MathUtils;

/**
 * The task of dividing up the reference database into partitions for
 * parallel processing.
 * @author mullermw
 *
 */
public class ReferencePartitioningTask implements Task, Serializable {

	public static final long serialVersionUID = 1;
	
	private static Logger logger = Logger.getLogger(ReferencePartitioningTask.class.getSimpleName());
	static {
		logger.setLevel(Level.FINE);
	}
	
	/*inputs */
	private File referenceFile;
	private long maxMemory;
	private File dest;
	private int readLength;
	private File exonReference;
	private transient JobSubmissionParameters jobSubmissionParameters;
	
	/*outputs */
	private transient Boolean done = null;
	private transient int numberOfPartitions = -1;
	private transient List<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = null;
	
	private transient boolean remote;
	
	public ReferencePartitioningTask(File referenceFile, long maxMemory, File dest, int readLength, File exonReference) {
		if (referenceFile == null) throw new IllegalArgumentException("referenceFile cannot be null.");
		if (maxMemory < 1000000) throw new IllegalArgumentException("maxMemory cannot be less than 1000000");
		if (dest == null) throw new IllegalArgumentException("dest cannot be null.");
		if (!dest.isDirectory()) throw new IllegalArgumentException("dest is not a directory.");
		if (!dest.canWrite()) throw new IllegalArgumentException("dest is not writeable.");
		if (readLength < 25) throw new IllegalArgumentException("readLength cannot be less than 25");
		this.referenceFile = referenceFile == null ? null : referenceFile.getAbsoluteFile();
		this.maxMemory = maxMemory;
		this.dest = dest == null ? null : dest.getAbsoluteFile();
		this.readLength = readLength;
		this.exonReference = exonReference == null ? null : exonReference.getAbsoluteFile();
	}
	
	@Override
	public void doTask() throws TaskException {
    	FileUtils.deleteDir(getWorkingDir());
    	File[] filesToDelete = dest.listFiles(new FileFilter() {
    		@Override
    		public boolean accept(File pathname) {
    			if (pathname.getName().equals(getJunctionFastaFile().getName())) return true;
    			if (pathname.getName().equals(getJunctionInfoFile().getName())) return true;
    			if (pathname.getName().startsWith("reference") && !pathname.getName().endsWith(".idx")) return true;
    			return false;
    		}
    	});
    	for (File file : filesToDelete)
    		file.delete();
		if (this.remote)
			doTaskRemotely();
		else
			doTaskLocally();
	}
	
	@Override
	public boolean isRemote() {
		return remote;
	}
	
	@Override
	public void setRemote(boolean remote) {
		this.remote = remote;
	}
	
	public static void main(String[] args) throws IOException, ClassNotFoundException, TaskException {
		File paramFile = new File(args[0]);
		ObjectInputStream stream = null;
		try {
			stream = new ObjectInputStream(new FileInputStream(paramFile));
			ReferencePartitioningTask task = (ReferencePartitioningTask)stream.readObject();
			task.doTaskLocally();
		} finally {
			if (stream != null) stream.close();
		}
	}
	
	public void doTaskRemotely() throws TaskException {
		ObjectOutputStream stream = null;
		try {
			File paramFile = getParamFile();
			paramFile.createNewFile();
			stream = new ObjectOutputStream(new FileOutputStream(paramFile));
			stream.writeObject(this);
			stream.close();
			String[] cmdArr = {
				"java", "-Xmx1g", "-classpath", getClassPath(), this.getClass().getName(), getParamFile().getAbsolutePath()
			};
			String cmd = Utilities.join(" ", cmdArr);
			ArrayList<String> cmds = new ArrayList<String>();
			cmds.add(cmd);
			jobSubmissionParameters.setMemoryRequirement(2 * MathUtils.GIGABYTE);
			ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(jobSubmissionParameters);
			clusterInterface.executeJob(getClusterScript(), cmds);
	        logger.info("Waiting for reference partitioning jobs to finish...");
	        while (!clusterInterface.checkIfLoggedJobsComplete())
	            Thread.sleep(5000);
	    	
	        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
	        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
	        	StringBuffer msg = new StringBuffer("The following jobs failed:\n");
	        	for (File file : listOfScriptOutputFilesForFailedJobs)
	                msg.append(file.getPath()+"\n");
	        	logger.severe(msg.toString());
	            throw new Exception("One or more reference partitioning jobs failed.");
	        }
	        
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
			throw new TaskException(e);
		} finally {
			try {
				if (stream != null) stream.close();
			} catch (IOException e) {
				throw new TaskException(e);
			}
		}
		if (!isDone()) throw new TaskException("Ref partitioning has failed.");
	}
	
	public void doTaskLocally() throws TaskException {
		this.numberOfPartitions = -1;
		done = false;
		PrintWriter tableWriter = null;
        try {
        	IndexedFastaDatabase fastaDBOfFullReference = new IndexedFastaDatabase(referenceFile, new File(dest, "reference.idx"));
        	File junctionFastaFile = getJunctionFastaFile();

			if (exonReference != null) {
				logger.info("Extracting splice junctions");
				//Extract Splice Junctions
				InputStream exonReferenceStream = null;
				PrintStream junctionStream = null;
				PrintStream infStream = null;
				try {
					exonReferenceStream = new FileInputStream(exonReference);
					File infFile = getJunctionInfoFile();
					junctionStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(junctionFastaFile, false/*append*/)), true/*autoflush*/);
					
					SpliceJunctionExtractor extractor = new SpliceJunctionExtractor(fastaDBOfFullReference, exonReferenceStream, junctionStream, getFlankSize());
					long junctionCount = extractor.extractJunctions();
					junctionStream.close();
					infStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(infFile, false)), true);
					Properties header = new Properties();
					header.put(FileHeader.CREATED_BY_KEY, SpliceJunctionExtractor.class.getSimpleName());
					header.put(FileHeader.SEQUENCE_SOURCE_FILE, referenceFile.getAbsolutePath());
					header.put(FileHeader.GENE_MODEL_SOURCE_FILE, exonReference.getAbsolutePath());
					header.put(FileHeader.JUNCTION_FLANK_SIZE, Integer.toString(getFlankSize()));
					header.put(FileHeader.JUNCTION_FASTA_FILE, junctionFastaFile.getPath());
					header.put(FileHeader.FILE_SIZE, Long.toString(junctionFastaFile.length()));
					header.put(FileHeader.JUNCTION_COUNT, Long.toString(junctionCount));
					FileHeader.writeHeader(header, infStream, "This file describes a file containing sequence flanking known and putative splice junctions.");

					infStream.close();
				} finally {
					if (junctionStream != null) junctionStream.close();
					if (exonReferenceStream != null) exonReferenceStream.close();
					if (infStream != null) infStream.close();
				}
			} else {
				logger.info("No exon reference provided.  Junctions will not be extracted.");
			}
				
			long maxRefSeqLen = ReadMapper.calculateMaxReferenceSizeForMapReads(maxMemory);
			File fileReferencePartitionTable = getReferencePartitionTableFile();
			
			logger.info("Started reference partitioning");
	        logger.info("Reference file will be partitioned into files with no more than " + (new DecimalFormat("0.00E0")).format(maxRefSeqLen) + " bases each.");
	        

        	tableWriter = new PrintWriter(new BufferedWriter(new FileWriter(fileReferencePartitionTable)));
	        @SuppressWarnings("unchecked")
	        ArrayList<HashSet<String>> listOfReferenceHeaderSets = fastaDBOfFullReference.getSetsOfSequenceHeadersEachTotalingLessThanNMegabases(maxRefSeqLen);
	        this.numberOfPartitions = listOfReferenceHeaderSets.size();
	        File[] filesReferenceFastaPartitions = new File[this.numberOfPartitions];
	        for (int i = 0; i < this.numberOfPartitions; i++)
	            filesReferenceFastaPartitions[i] = getReferencePartitionFile(i);
	        
	        Properties p = new Properties();
	        p.put(FileHeader.CREATED_BY_KEY, "FastaDatabase.partitionFastaFileToFastaFilesByHeaderSets()");
	        p.put(FileHeader.SEQUENCE_SOURCE_FILE, fastaDBOfFullReference.getFileFastaDatabase().getAbsolutePath());
	        p.put(FileHeader.MAX_MEMORY_KEY, Long.toString(this.maxMemory));
	        FileHeader.writeHeader(p, tableWriter, "Maps the entries in the reference partitions back to the original reference");
	        this.listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = FastaDatabase.partitionFastaFileToFastaFilesByHeaderSets(tableWriter, filesReferenceFastaPartitions, fastaDBOfFullReference, listOfReferenceHeaderSets);
        } catch (Exception e) {
        	logger.log(Level.SEVERE, "Failure in Reference Partitioning", e);
        	throw new TaskException("Failure in Reference Partitioning.", e);
        } finally {
        	if (tableWriter != null) tableWriter.close();
        }
        done = true;
        logger.info("Finished reference partitioning");

	}

	@Override
	public boolean isDone() throws TaskException {
		logger.info("Checking if reference partitioning has been successfully completed.");
		if (done != null) return done;
		
		this.numberOfPartitions = -1;
		this.listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = null;
		File tableFile = getReferencePartitionTableFile();
		logger.fine("looking for table file");
		if (!tableFile.exists()) return false;
		InputStream stream = null;
		try {
			try {
				this.listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = WholeTranscriptomeAnalyzer.loadPartitionTable(tableFile);
				stream = new FileInputStream(tableFile);
				Map<String, FileHeader.Value> header = FileHeader.parseHeader(stream);
				
				//Check the Partition Table.
				logger.fine("Comparing Source Files");
				FileHeader.Value value = header.get(FileHeader.SEQUENCE_SOURCE_FILE);
				if (value == null) return false;
				File sourceFile = new File(value.toString()).getAbsoluteFile();
				logger.finer(sourceFile.toString());
				logger.finer(this.referenceFile.getAbsoluteFile().toString());
				if (!this.referenceFile.getAbsoluteFile().equals(sourceFile)) return false;
				
				logger.fine("Comparing maxMemory");
				value = header.get(FileHeader.MAX_MEMORY_KEY);
				if (value == null || value.longValue() != this.maxMemory) return false;
				
				logger.fine("Checking for partition files");
				stream = new FileInputStream(tableFile);
				BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
				int numberOfPartitions = 0;
				for (String line = reader.readLine(); line != null; line = reader.readLine()) {
					if (line.matches("(\\d+)\\t\\d+\\t\\d+\\t.*")) {
						int index = Integer.parseInt(line.substring(0, line.indexOf("\t")));
						if (index + 1 > numberOfPartitions) numberOfPartitions = index + 1;
						File partitionFile = getReferencePartitionFile(index);
						logger.finer("Looking for "+partitionFile+".");
						if (!partitionFile.exists()) return false;
					}
				}
				this.numberOfPartitions = numberOfPartitions;
				//Check the junction Reference.
				logger.fine("Checking junction reference");
				if (exonReference != null) {
					File junctionInfoFile = getJunctionInfoFile();
					if (junctionInfoFile.exists() == false) return false;
					try {
						stream = new FileInputStream(junctionInfoFile);
						header = FileHeader.parseHeader(stream);
					} finally {
						if (stream != null) stream.close();
					}
					value = header.get(FileHeader.SEQUENCE_SOURCE_FILE);
					if (value == null) return false;
					File fileFromHeader = new File(value.toString()).getAbsoluteFile();
					if (!this.referenceFile.getAbsoluteFile().equals(fileFromHeader)) return false;
					
					value = header.get(FileHeader.GENE_MODEL_SOURCE_FILE);
					if (value == null) return false;
					fileFromHeader = new File(value.toString()).getAbsoluteFile();
					if (!this.exonReference.getAbsoluteFile().equals(fileFromHeader)) return false;
					
					value = header.get(FileHeader.JUNCTION_FLANK_SIZE);
					if (value == null) return false;
					if (this.getFlankSize() != value.intValue()) return false;
					
					File junctionFastaFile = new File(header.get(FileHeader.JUNCTION_FASTA_FILE).toString());
					if (junctionFastaFile.equals(this.getJunctionFastaFile()) == false) {
						logger.severe(String.format("Path Mismatch between %s and %s.", junctionFastaFile, junctionInfoFile));
						return false;
					}
					long fileSize = header.get(FileHeader.FILE_SIZE).longValue();
					if (junctionFastaFile.length() != fileSize) {
						logger.severe(String.format("File size mismatch between %s and %s.", junctionFastaFile, junctionInfoFile));
						return false;
					}
				}
				
			} finally {
				if (stream != null) stream.close();
			}
		} catch (Exception e) {
			logger.log(Level.SEVERE, "Failure in determining if Reference Partitioning is done.", e);
			throw new TaskException("Failure in determining if Reference Partitioning is done.", e);
		}
		logger.info("Reference partitioning looks ok.");
		return true;
	}

	public File getReferencePartitionTableFile() {
		return new File(dest, "reference.partition_table.out");
	}
	
	public File getReferencePartitionFile(int index) {
		return new File(dest, "/reference." + index + ".fasta");
	}
	
	public File getWorkingDir() {
		File f = new File(dest, "ref_partitioning");
		if (!f.exists()) f.mkdir();
		return f;
	}
	
	public File getClusterScript() {
		return new File(getWorkingDir(), "script.sh");
	}
	
	public File[] getPartitions() {
		File[] arr = new File[getNumberOfPartitions()];
		for (int i=0; i<arr.length; i++) {
			arr[i] = getReferencePartitionFile(i);
		}
		return arr;
	}
	
	public int getNumberOfPartitions() {
		return this.numberOfPartitions;
	}

	public File getJunctionFastaFile() {
		return new File(dest, "junctions.fa");
	}
	
	public File getJunctionInfoFile() {
		return new File(dest, "junctions.inf");
	}
		
	public int getFlankSize() {
		return readLength - 4;
	}

	public List<Integer>[] getListsOfFullRefSequenceNumbersOrderedForEachReferencePartition() {
		return listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
	}

	public void setListsOfFullRefSequenceNumbersOrderedForEachReferencePartition(
			List<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition) {
		this.listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
	}
	
	public JobSubmissionParameters getJobSubmissionParameters() {
		return jobSubmissionParameters;
	}
	public void setJobSubmissionParameters(JobSubmissionParameters jobSubmissionParameters) {
		this.jobSubmissionParameters = jobSubmissionParameters;
	}
	public String getClassPath() {
		return WholeTranscriptomeAnalyzer.getLibDir().getPath() + "/'*':"+WholeTranscriptomeAnalyzer.getWholeTranscriptomeJar().getPath();
	}
	public File getParamFile() {
		return new File(this.getWorkingDir(), "params.obj");
	}
}

