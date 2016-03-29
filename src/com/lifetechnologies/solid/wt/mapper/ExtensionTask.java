package com.lifetechnologies.solid.wt.mapper;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.CompressionUtilities;
import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.ReadMapper;
import com.lifetechnologies.solid.wt.ReadMappingExtender;
import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.TextFileUtilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.util.FileUtils;

/**
 * The task of extending the mapreads anchors.
 * @author mullermw
 *
 */
public class ExtensionTask implements Task {
	
	private Logger logger = Logger.getLogger(ExtensionTask.class.getName());
	
	private long maxBytesMemoryPerJob;
	private double memoryRequirementAdjustmentFactor;
	private JobSubmissionParameters parameterTemplate;
	private File fileFullLengthReadsCsfasta;
	private int positionOfLeftRead;
	private File[] maFilesLeft;
	private int positionOfRightRead;
	private File[] maFilesRight;
	private File maFileJunctions;
	private File[] filesReferenceFastaPartitions;
	private File junctionReference;
	private ReadMask readMask;
	private Boolean applyValidAdjacentRules;
	private Boolean matchIubCodes;
	private int matchScore = 1;
	private int mismatchPenalty = -1;
	private File outputDir;
	private File scratchDir;
	private Boolean deleteIntermediateFiles;
	private Boolean compressIntermediateFiles;
	
	private int numberOfReadsFilePartitions = 1;
	
	@Override
	public void doTask() throws TaskException {

        logger.info("Started extension jobs");
        try {
        	
        	FileUtils.deleteFiles(getTmpDir(), new FilenameFilter() {
        		@Override
        		public boolean accept(File dir, String name) {
        			if (name.endsWith(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE)) return true;
        			if (name.endsWith(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE+WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE)) return true;
        			if (name.endsWith(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE+".sh")) return true;
        			if (name.endsWith(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE+".sh.out" )) return true;
        			return false;
        		}
        	});
        	getMasterJobRemovalScript().delete();
        	getMasterJobSubmissionScript().delete();
        	getMasterListOfScriptOutputFiles().delete();
        	
	        ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();
	
	        long memoryRequiredByExtension = (long)Math.ceil(ReadMappingExtender.calculateMemoryRequiredByExtendMapReadsInBytes(ReadMapper.calculateMaxReferenceSizeForMapReads(maxBytesMemoryPerJob))*memoryRequirementAdjustmentFactor);
	        
	        JobSubmissionParameters params = parameterTemplate.clone();
	        params.setMemoryRequirement(memoryRequiredByExtension);
	        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
	
	        int readLength = new FastaDatabase(fileFullLengthReadsCsfasta).getMonoLength(10).intValue() - 1;
	        
	        for (int indexPartition = 0; indexPartition<filesReferenceFastaPartitions.length; indexPartition++) {
	        	File referencePartition = filesReferenceFastaPartitions[indexPartition];
	        	File maFileLeft = maFilesLeft[indexPartition];
	        	File maFileRight = null;
	        	if (maFilesRight.length > 0) maFileRight= maFilesRight[indexPartition];
	        	File maxFileLeft = getMaxForMaFile(maFileLeft);
	        	File maxFileRight = getMaxForMaFile(maFileRight);
	        	File mergedMaxFileLeft = getMergedMaxForMaFile(maFileLeft);
	        	File mergedMaxFileRight = getMergedMaxForMaFile(maFileRight);
	        	
	        	for (File file : new File[] { maFileLeft, maFileRight, maFileJunctions }) {
	        		if (file != null && !file.exists()) {
	        			File compressedFile = new File(file.getParent(), file.getName()+WholeTranscriptomeAnalyzer.EXTENSION_FOR_COMPRESSED_FILE);
	        			if (compressedFile.exists()) {
	        				logger.info("Decompressing intermediate files");
	        				CompressionUtilities.decompressFiles(compressedFile, compressedFile.getParentFile());
	        				logger.info("OK");
	        			} else {
	                        throw new Exception("Could not find MA file:\t" + file.getPath() + "\nPlease re-run with RUN_MAPPING turned on.");
	        			}
	        		}
	        	}
	        	
	        	for (File file : new File[] { mergedMaxFileLeft, mergedMaxFileRight} )
	        		if (file != null && file.exists()) file.delete();
	        	
	        	ReadMappingExtender extender = new ReadMappingExtender(WholeTranscriptomeAnalyzer.getExtendMappedReadsPy(), referencePartition,
	                    maFileLeft, fileFullLengthReadsCsfasta, readLength);

	        	extender.startExtendMappedReadsJob(matchScore, mismatchPenalty,
	                    this.applyValidAdjacentRules,
	                    this.matchIubCodes,
	                    positionOfLeftRead,
	                    clusterInterface, "wt_extension.",
	                    this.scratchDir, maFileLeft.getParentFile(),
	                    maxFileLeft, new File(maFileLeft.getParent(), maxFileLeft.getName()+WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE),
	                    true, this.readMask.asBitString(readLength));
	        	
	        	listOfFilesThatAreNoLongerNeeded.add(maFileLeft);
	        	
	        	if (maFileRight != null) {
		        	extender = new ReadMappingExtender(WholeTranscriptomeAnalyzer.getExtendMappedReadsPy(), referencePartition,
		                    maFileRight, fileFullLengthReadsCsfasta, readLength);
		        	
		        	extender.startExtendMappedReadsJob(matchScore, mismatchPenalty,
		                    this.applyValidAdjacentRules,
		                    this.matchIubCodes,
		                    positionOfRightRead,
		                    clusterInterface, "wt_extension.",
		                    this.scratchDir, maFileRight.getParentFile(),
		                    maxFileRight, new File(maFileRight.getParent(), maxFileRight.getName()+WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE),
		                    true, this.readMask.asBitString(readLength));
		        	
		        	listOfFilesThatAreNoLongerNeeded.add(maFileRight);
	        	}
	        }
	        
	        if (junctionReference != null) {
	        	File junctionMaxFile = getJunctionMaxFile();
	        	File indexFile = new File(junctionMaxFile.getParent(), junctionMaxFile.getName()+WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE);
	        	
	        	ReadMappingExtender extender = new ReadMappingExtender(WholeTranscriptomeAnalyzer.getExtendMappedReadsPy(), junctionReference, maFileJunctions, fileFullLengthReadsCsfasta, readLength);
	        	extender.startExtendMappedReadsJob(this.matchScore, this.mismatchPenalty, this.applyValidAdjacentRules, this.matchIubCodes, positionOfLeftRead, clusterInterface, "wt_extension.", this.scratchDir, this.maFileJunctions.getParentFile(), junctionMaxFile, indexFile, true, this.readMask.asBitString(readLength));
	        	
	        	listOfFilesThatAreNoLongerNeeded.add(maFileJunctions);
	        }
	        
	        
	        File masterSubmissionScript = getMasterJobSubmissionScript();
	        clusterInterface.writeMasterJobSubmissionFileFromLog(masterSubmissionScript);
	        File masterRemovalScript = getMasterJobRemovalScript();
	        clusterInterface.writeMasterJobRemovalFileFromLog(masterRemovalScript);
	        clusterInterface.writeMasterListOfJobOutputFilesFromLog(getMasterListOfScriptOutputFiles());
	        
	        masterSubmissionScript.setExecutable(true);
	        masterRemovalScript.setExecutable(true);
	        
	     // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
	        logger.info("Waiting for extension jobs to finish...");
	        while (!clusterInterface.checkIfLoggedJobsComplete())
	            Thread.sleep(5000);
	
	        logger.info("OK");
	        
	        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
	        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
	        	StringBuffer msg = new StringBuffer("The following jobs failed:\n");
	            Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
	            while (iteratorOverFailFiles.hasNext())
	                msg.append(( iteratorOverFailFiles.next()).getPath()+"\n");
	            logger.severe(msg.toString());
	            throw new TaskException("One or more extension jobs failed.");
	        }
	        
	        if (deleteIntermediateFiles || compressIntermediateFiles) {
	            if (compressIntermediateFiles && deleteIntermediateFiles)  logger.info("Compressing and deleting intermediate files...");
	            else if (compressIntermediateFiles)   logger.info("Compressing intermediate files...");
	            else    logger.info("Deleting intermediate files...");
	            for (File file : listOfFilesThatAreNoLongerNeeded) {
	                if (compressIntermediateFiles)
	               		CompressionUtilities.compressFile(file, new File(file.getPath() + WholeTranscriptomeAnalyzer.EXTENSION_FOR_COMPRESSED_FILE));
	                if (deleteIntermediateFiles)
	                    file.delete();
	            }
	            logger.info("OK");
	        }
	
	
	        logger.info("Finished extension jobs ");
        } catch (Exception e) {
        	logger.log(Level.SEVERE, e.getMessage(), e);
        	if (e instanceof TaskException) throw (TaskException)e;
        	throw new TaskException(e.getMessage(), e);
        }
	}
	
	@Override
	public boolean isDone() throws TaskException {
        logger.info("Checking that extension jobs from previous run completed successfully...");
        try {
        	WholeTranscriptomeAnalyzer.checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(getMasterListOfScriptOutputFiles(), "\t", 0));
        } catch (Exception e) {
        	return false;
        }
        logger.info(" OK");
		return true;
	}
	
	@Override
	public boolean isRemote() {
		return true;
	}
	
	@Override
	public void setRemote(boolean remote) {
		if (remote == false) throw new IllegalArgumentException("ExtensionTask only supported as a remote task.");
	}
	
	private static File getMaxForMaFile(File maFile) {
		if (maFile == null ) return null;
		return new File(maFile.getParent(), maFile.getName().replaceFirst(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE +"$", WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE));
	}
	
	private static File getMergedMaxForMaFile(File maFile) {
		if (maFile == null) return null;
		return new File(maFile.getParent(), maFile.getName().replaceFirst(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE+"$", WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE));
	}
	
	public File getScriptDir() {
		return new File(getOutputDir(), "scripts");
	}
	
	public File getTmpDir() {
		return new File(getOutputDir(), "tmp");
	}
	
	public File getMasterJobSubmissionScript() {
		return new File (getScriptDir(), "extension.master.script_submission.sh");
	}
	
	public File getMasterJobRemovalScript() {
		return new File(getScriptDir(), "extension.master.script_removal.sh");
	}
	
	public File getMasterListOfScriptOutputFiles() {
		return new File(getOutputDir(), "tmp/extension.list_of_script_output_files.out");
	}
	
	public File getJunctionMaxFile() {
		if (this.maFileJunctions == null) return null;
		return new File(maFileJunctions.getParent(), maFileJunctions.getName().replaceAll(WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE+"$", WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE));
	}

	public long getMaxBytesMemoryPerJob() {
		return maxBytesMemoryPerJob;
	}

	public void setMaxBytesMemoryPerJob(long maxBytesMemoryPerJob) {
		this.maxBytesMemoryPerJob = maxBytesMemoryPerJob;
	}

	public double getMemoryRequirementAdjustmentFactor() {
		return memoryRequirementAdjustmentFactor;
	}

	public void setMemoryRequirementAdjustmentFactor(
			double memoryRequirementAdjustmentFactor) {
		this.memoryRequirementAdjustmentFactor = memoryRequirementAdjustmentFactor;
	}

	public JobSubmissionParameters getParameterTemplate() {
		return parameterTemplate;
	}

	public void setParameterTemplate(JobSubmissionParameters parameterTemplate) {
		this.parameterTemplate = parameterTemplate;
	}

	public File getFileFullLengthReadsCsfasta() {
		return fileFullLengthReadsCsfasta;
	}

	public void setFileFullLengthReadsCsfasta(File fileFullLengthReadsCsfasta) {
		this.fileFullLengthReadsCsfasta = fileFullLengthReadsCsfasta;
	}

	public int getPositionOfLeftRead() {
		return positionOfLeftRead;
	}

	public void setPositionOfLeftRead(int positionOfLeftRead) {
		this.positionOfLeftRead = positionOfLeftRead;
	}

	public File[] getMaFilesLeft() {
		return maFilesLeft;
	}

	public void setMaFilesLeft(File[] maFilesLeft) {
		this.maFilesLeft = maFilesLeft;
	}

	public int getPositionOfRightRead() {
		return positionOfRightRead;
	}

	public void setPositionOfRightRead(int positionOfRightRead) {
		this.positionOfRightRead = positionOfRightRead;
	}

	public File[] getMaFilesRight() {
		return maFilesRight;
	}

	public void setMaFilesRight(File[] maFilesRight) {
		if (maFilesRight == null) maFilesRight = new File[0];
		this.maFilesRight = maFilesRight;
	}

	public File getMaFileJunctions() {
		return maFileJunctions;
	}

	public void setMaFileJunctions(File maFileJunctions) {
		this.maFileJunctions = maFileJunctions;
	}

	public File[] getFilesReferenceFastaPartitions() {
		return filesReferenceFastaPartitions;
	}

	public void setFilesReferenceFastaPartitions(
			File[] filesReferenceFastaPartitions) {
		this.filesReferenceFastaPartitions = filesReferenceFastaPartitions;
	}

	public File getJunctionReference() {
		return junctionReference;
	}

	public void setJunctionReference(File junctionReference) {
		this.junctionReference = junctionReference;
	}

	public ReadMask getReadMask() {
		return readMask;
	}

	public void setReadMask(ReadMask readMask) {
		this.readMask = readMask;
	}

	public Boolean getApplyValidAdjacentRules() {
		return applyValidAdjacentRules;
	}

	public void setApplyValidAdjacentRules(Boolean applyValidAdjacentRules) {
		this.applyValidAdjacentRules = applyValidAdjacentRules;
	}

	public Boolean getMatchIubCodes() {
		return matchIubCodes;
	}

	public void setMatchIubCodes(Boolean matchIubCodes) {
		this.matchIubCodes = matchIubCodes;
	}
	
	public int getMatchScore() {
		return matchScore;
	}

//	public void setMatchScore(int matchScore) {
//		this.matchScore = matchScore;
//	}

	public int getMismatchPenalty() {
		return mismatchPenalty;
	}

//	public void setMismatchPenalty(int mismatchPenalty) {
//		this.mismatchPenalty = mismatchPenalty;
//	}

	public File getOutputDir() {
		return outputDir;
	}

	public void setOutputDir(File outputDir) {
		this.outputDir = outputDir;
	}

	public File getScratchDir() {
		return scratchDir;
	}

	public void setScratchDir(File scratchDir) {
		this.scratchDir = scratchDir;
	}

	public Boolean getDeleteIntermediateFiles() {
		return deleteIntermediateFiles;
	}

	public void setDeleteIntermediateFiles(Boolean deleteIntermediateFiles) {
		this.deleteIntermediateFiles = deleteIntermediateFiles;
	}

	public Boolean getCompressIntermediateFiles() {
		return compressIntermediateFiles;
	}

	public void setCompressIntermediateFiles(Boolean compressIntermediateFiles) {
		this.compressIntermediateFiles = compressIntermediateFiles;
	}

	public int getNumberOfReadsFilePartitions() {
		return numberOfReadsFilePartitions;
	}

	public void setNumberOfReadsFilePartitions(int numberOfReadsFilePartitions) {
		throw new UnsupportedOperationException("setNumberOfReadsFilePartitions is not supported.");
	}
}
