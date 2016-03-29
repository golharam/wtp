package com.lifetechnologies.solid.wt.mapper;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.util.FileHeader;
import com.lifetechnologies.util.MathUtils;

/**
 * The task of 'splitting' the reads into pieces for anchoring.
 * @author mullermw
 *
 */
public class SplitReadsTask implements Task, Serializable {
	
	public static final long serialVersionUID = 1;
	
	/* Input params */
	private File readsFile;
	private File destDir;
	private int readLength = -1;
	private int lengthLeft;
	private int lengthRight;
	
	
	/* Output params */
	private transient int leftReadStartsAt = 1;
	private transient int rightReadStartsAt;
	private transient int readCount = -1;
	private transient Boolean done = null;
	
	private transient JobSubmissionParameters jobSubmissionParameters;
	private transient boolean remote;
	
	private static Logger logger = Logger.getLogger(SplitReadsTask.class.toString());
	
	public SplitReadsTask(File readsFile, File destDir, int lengthLeft, int lengthRight) throws ParseException {
		if (lengthLeft < 15) throw new IllegalArgumentException("lengthLeft cannot be less than 15");
		if (lengthRight < 0) throw new IllegalArgumentException("lengthRight cannot be less than 0");
		if (readsFile == null) throw new IllegalArgumentException("reads cannot be null");
		if (!readsFile.exists()) throw new IllegalArgumentException("readsFile does not exist.");
		if (!readsFile.canRead()) throw new IllegalArgumentException("readsFile is not readable.");
		if (destDir == null) throw new IllegalArgumentException("destDir cannot be null");
		if (!destDir.isDirectory()) throw new IllegalArgumentException("destDir is not a directory.");
		if (!destDir.canRead()) throw new IllegalArgumentException("destDir is not readable.");
		if (!destDir.canWrite()) throw new IllegalArgumentException("destDir is not writeable.");
		
		Map<String, Long> readMap = null;
		try {
			FastaDatabase db = new FastaDatabase(readsFile);
			readMap = db.getSequenceLengths(20);
		} catch (Exception e) {
			throw new ParseException(e.getMessage(), 1);
		}
		for (Long length : readMap.values()) {
			if (length > 200) throw new IllegalArgumentException("Reads are too long.");
			if (this.readLength < 0) this.readLength = length.intValue() - 1;
			else if (this.readLength != length.intValue() - 1) throw new IllegalArgumentException("This reads stream contains reads of varying length.");
		}
		
		if (lengthLeft > readLength) throw new IllegalArgumentException("lengthLeft cannot be greater than readLength.");
		if (lengthRight > readLength - 1 ) throw new IllegalArgumentException("lengthRight must be less than readLength - 1.");
		
		this.readsFile = readsFile.getAbsoluteFile();
		this.destDir = destDir.getAbsoluteFile();
		this.lengthLeft = lengthLeft;
		this.lengthRight = lengthRight;
	}
	
	@Override
	public void doTask() throws TaskException {
        File leftReadsFile = getLeftReadsFile();
        File rightReadsFile = getRightReadsFile();
        if (leftReadsFile.exists()) leftReadsFile.delete();
        if (rightReadsFile != null && rightReadsFile.exists()) rightReadsFile.delete();
        if (this.isRemote()) {
        	doTaskRemotely();
        } else {
        	doTaskLocally();
        }
	}
	
	public static void main(String[] args) throws IOException, ClassNotFoundException, TaskException {
		File paramFile = new File(args[0]);
		ObjectInputStream stream = null;
		try {
			stream = new ObjectInputStream(new FileInputStream(paramFile));
			SplitReadsTask task = (SplitReadsTask)stream.readObject();
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
	        logger.info("Waiting for splitting jobs to finish...");
	        while (!clusterInterface.checkIfLoggedJobsComplete())
	            Thread.sleep(5000);
	    	
	        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
	        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
	        	StringBuffer msg = new StringBuffer("The following jobs failed:\n");
	        	for (File file : listOfScriptOutputFilesForFailedJobs)
	                msg.append(file.getPath()+"\n");
	        	logger.severe(msg.toString());
	            throw new Exception("One or more splitting jobs failed.");
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
		done = false;
		logger.info("Started read splitting");
        logger.info("\nSplitting the read sequences:\n" +
                " -->First part of sequence will be colors 1 to " + lengthLeft + ":\n" +
                " -->Last part of sequence will be colors " + (readLength - lengthRight +1) + " to " + readLength );
        try {
        	FastaDatabase reads = new FastaDatabase(readsFile);
        	reads.generateSplitVersionsToFiles(lengthLeft +1, lengthRight, "T0", getLeftReadsFile(), getRightReadsFile());
        	this.readCount = reads.getNumberOfSequencesInDatabase();
        } catch (Exception e) {
        	throw new TaskException(e.getMessage(), e);
        }
        
        logger.info("Finished read splitting");

		// the position of the last part of a read is shifted down by one from the
		// true position because we must add the faux "0" color at the front
		rightReadStartsAt = readLength - lengthRight;
		done = true;
	}
	
	public boolean isDone() throws TaskException {
		logger.info("Checking if Splitting has successfully completed");
		if (done != null) return done;
		File leftReadsFile = getLeftReadsFile();
		File rightReadsFile = getRightReadsFile();
		if (!leftReadsFile.exists()) {
			logger.warning("No Left Reads File.");
			return false;
		}
		if (rightReadsFile != null && !rightReadsFile.exists()) {
			logger.warning("No Right Reads File.");
			return false;
		}
		InputStream leftSplitStream = null;
		InputStream rightSplitStream = null;
		try {
			leftSplitStream = new FileInputStream(leftReadsFile);
			Map<String, FileHeader.Value> props = FileHeader.parseHeader(leftSplitStream);
			
			FileHeader.Value value = props.get(FileHeader.SEQUENCE_SOURCE_FILE);
			if (value == null) {
				logger.warning("No Sequence Source File specified in header of "+leftReadsFile );
				return false;
			}
			File sourceFile = new File(value.toString());
			if (!sourceFile.getAbsoluteFile().equals(this.readsFile)) {
				logger.warning("File path discrepancy.");
				return false;
			}
			
			value = props.get(FileHeader.SPLIT_KEY);
			if (value == null || !value.toString().trim().toLowerCase().equals("left")) {
				logger.warning("incorrect header in "+leftReadsFile);
				return false;
			}
			
			value = props.get(FileHeader.SPLIT_LENGTH_KEY);
			if (value == null || !value.isNumeric() || !value.intValue().equals(lengthLeft)) {
				logger.warning("Read length mismatch in header");
				return false;
			}
			
			if (rightReadsFile == null) return true;
			
			rightSplitStream = new FileInputStream(rightReadsFile);
			props = FileHeader.parseHeader(rightSplitStream);
			
			value = props.get(FileHeader.SEQUENCE_SOURCE_FILE);
			if (value == null)  {
				logger.warning("No Source File specified in header.");
				return false;
			}
			sourceFile = new File(value.toString());
			if (!sourceFile.getAbsoluteFile().equals(this.readsFile)) {
				logger.warning("File path discrepancy.");
				return false;
			}
			
			value = props.get(FileHeader.SPLIT_KEY);
			if (value == null || !value.toString().trim().toLowerCase().equals("right")) { 
				logger.warning("incorrect header in " + rightReadsFile);
				return false;
			}
			
			value = props.get(FileHeader.SPLIT_LENGTH_KEY);
			if (value == null || !value.isNumeric() || !value.intValue().equals(lengthRight)) {
				logger.warning("Read length mismatch in header.");
				return false;
			}
			
		} catch (IOException e) {
			throw new TaskException(e);
		} finally {
			try {
				if (leftSplitStream != null) leftSplitStream.close();
				if (rightSplitStream != null) rightSplitStream.close();
			} catch (IOException e) {
				throw new TaskException(e);
			}
		}
		return true;
	} 
	
	@Override
	public boolean isRemote() {
		return this.remote;
	}
	
	@Override
	public void setRemote(boolean remote) {
		this.remote = remote;
	}

	public File getLeftReadsFile() {
		return new File(destDir, "reads.first_"+lengthLeft+".csfasta");
	}

	public File getRightReadsFile() {
		if (lengthRight < 1) return null;
		return new File(destDir, "reads.last_"+lengthRight+".csfasta");
	}
	
	public File[] getSplitReadsFiles() {
		File leftFile = getLeftReadsFile();
		File rightFile = getRightReadsFile();
		if (rightFile == null)
			return new File[] { leftFile };
		else
			return new File[] { leftFile, rightFile};
	}
	
	public File getWorkingDir() {
		File dir = new File(destDir, "splitting");
		if (!dir.exists()) dir.mkdir();
		return dir;
	}
	
	public File getParamFile() {
		return new File(getWorkingDir(), "params.obj");
	}
	
	public File getClusterScript() {
		return new File(getWorkingDir(), "script.sh");
	}

	public String getClassPath() {
		return WholeTranscriptomeAnalyzer.getLibDir().getPath() + "/'*':"+WholeTranscriptomeAnalyzer.getWholeTranscriptomeJar().getPath();
	}
	
	public File getReadsFile() {
		return readsFile;
	}
	
	public int getReadLength() {
		return readLength;
	}

	public int getLengthLeft() {
		return lengthLeft;
	}

	public int getLengthRight() {
		return lengthRight;
	}

	public int getLeftReadStartsAt() {
		return leftReadStartsAt;
	}

	public int getRightReadStartsAt() {
		return rightReadStartsAt;
	}

	public int getReadCount() {
		return readCount;
	}

	public Logger getLogger() {
		return logger;
	}
	
	public JobSubmissionParameters getJobSubmissionParameters() {
		return jobSubmissionParameters;
	}

	public void setJobSubmissionParameters(
			JobSubmissionParameters jobSubmissionParameters) {
		this.jobSubmissionParameters = jobSubmissionParameters;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return this.getClass().getSimpleName()+":reads="+getReadsFile()+
		       ", leftReadsFile="+getLeftReadsFile()+
		       ", rightReads="+getRightReadsFile()+
		       ", readLength="+getReadLength() + 
		       ", lengthLeft="+getLengthLeft() +
		       ", lengthRight="+getLengthRight() +
		       ", leftReadStartsAt="+getLeftReadStartsAt()+
		       ", rightReadStartsAt="+getRightReadStartsAt()+
		       ", readCount="+getReadCount();
	}
}

