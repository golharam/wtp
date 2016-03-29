package com.lifetechnologies.solid.wt.mapper;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.ReadMapper;
import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.TextFileUtilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.util.FileUtils;

/**
 * The task of using mapreads to anchor the reads on the reference.
 * @author mullermw
 *
 */
public class MappingTask implements Task {

	private Logger logger = Logger.getLogger(MappingTask.class.getSimpleName());
	
	/* input */
	private long maxBytesMemoryPerJob;
	private double memoryRequirementAdjustmentFactor;
	private JobSubmissionParameters parameterTemplate;
	private File[] filesSplitReadsCSFasta;
	private ReadMask[] masksByReadSplit;
	private int[] allowedMismatches;
	private File[] filesReferenceFastaPartitions;
	private File junctionReference;
	private int numberOfReads;
	private Boolean applyValidAdjacentRules;
	private Boolean matchIubCodes;
	private int maxMatchingLocations;
	private File outputDir;
	private File scratchDir;
	
	private int numberOfReadsFilePartitions = 1;

	/* output */
	private File[] leftMaFiles;
	private File[] rightMaFiles;
	private File junctionMaFile;
	
	private int[] splitLengths;
	
	@Override
	public void doTask() throws TaskException {
		try {
			
        	FileUtils.deleteDir(getMappingAndExtensionFolder());
        	getMasterJobRemovalScript().delete();
        	getMasterJobSubmissionScript().delete();
        	getMasterListOfScriptOutputFiles().delete();
			
			checkParams();
			
			logger.info("Starting mapping jobs");
			long maxRefBasesPerJob = ReadMapper.calculateMaxReferenceSizeForMapReads(maxBytesMemoryPerJob);
			long memoryRequiredForMapReads = (long)Math.ceil(ReadMapper.calculateMemoryRequiredByMapReadsInBytes(maxRefBasesPerJob)*memoryRequirementAdjustmentFactor);
			JobSubmissionParameters params = parameterTemplate.clone();
			params.setMemoryRequirement(memoryRequiredForMapReads);
			ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
			leftMaFiles = new File[filesReferenceFastaPartitions.length];
			rightMaFiles = new File[filesReferenceFastaPartitions.length];
			for (int indexReadSplit = 0; indexReadSplit < filesSplitReadsCSFasta.length; indexReadSplit++) {
	
				for (int indexRefFile = 0; indexRefFile < filesReferenceFastaPartitions.length; indexRefFile++) {
	
					File subfolderForOutputFiles = new File(getMappingAndExtensionFolder() + "/reads_" + indexReadSplit + "/ref_" + indexRefFile);
					if (!subfolderForOutputFiles.exists())
						subfolderForOutputFiles.mkdirs();
					File splitReadsFile = filesSplitReadsCSFasta[indexReadSplit];
					int splitLength = this.splitLengths[indexReadSplit];
					ReadMapper mapper = new ReadMapper(WholeTranscriptomeAnalyzer.getMapreadsExe(), filesReferenceFastaPartitions[indexRefFile], splitReadsFile,
							splitLength, masksByReadSplit[indexReadSplit].asBitString(splitLength));
	
					int numberOfReadsPerJob = (int)Math.ceil((double)numberOfReads / numberOfReadsFilePartitions);
					for (int indexOfRead = 0; indexOfRead < numberOfReads; indexOfRead += numberOfReadsPerJob) {
	
						int maxReadIndexInCurrentRange = Math.min(indexOfRead + numberOfReadsPerJob, numberOfReads) -1;
	
						File fileMA = getMaFile(indexReadSplit, indexRefFile, splitLength);
						
						if (indexReadSplit == 0) {
							this.leftMaFiles[indexRefFile] = fileMA;
						} else if (indexReadSplit == 1) {
							this.rightMaFiles[indexRefFile] = fileMA;
						} 
						mapper.startMapReadsJob(true, indexOfRead, maxReadIndexInCurrentRange,
								allowedMismatches[indexReadSplit], applyValidAdjacentRules, matchIubCodes,
								this.maxMatchingLocations,
								clusterInterface, "wt_map.",
								this.scratchDir, WholeTranscriptomeAnalyzer.getSchemaDir(), subfolderForOutputFiles, fileMA, true);
					}
				}
			}

			if (junctionReference != null) {
				int splitLength = new FastaDatabase(filesSplitReadsCSFasta[0]).getMonoLength(10).intValue() - 1;
				ReadMapper mapper = new ReadMapper(WholeTranscriptomeAnalyzer.getMapreadsExe(), junctionReference, filesSplitReadsCSFasta[0],
												   splitLength, masksByReadSplit[0].asBitString(splitLength));
				File junctionOutputDir = getJunctionOutputDir();
				junctionOutputDir.mkdirs();
				if (!junctionOutputDir.exists()) throw new TaskException("Couldn't create directory:" + junctionOutputDir);

				this.junctionMaFile = getJunctionMaFile();
				mapper.startMapReadsJob(true, 0, this.numberOfReads - 1, 
						this.allowedMismatches[0], this.applyValidAdjacentRules, 
						this.matchIubCodes, this.maxMatchingLocations, 
						clusterInterface, "wt_map", this.scratchDir, 
						WholeTranscriptomeAnalyzer.getSchemaDir(), 
						junctionOutputDir, this.junctionMaFile, true);
			}
			
			File fileMasterMappingJobSubmissionScript = getMasterJobSubmissionScript();
			clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterMappingJobSubmissionScript);
			File fileMasterMappingJobRemovalScript = getMasterJobRemovalScript();
			clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterMappingJobRemovalScript);
			clusterInterface.writeMasterListOfJobOutputFilesFromLog(getMasterListOfScriptOutputFiles());
	
			fileMasterMappingJobSubmissionScript.setExecutable(true, true);
			fileMasterMappingJobRemovalScript.setExecutable(true, true);
	
			// FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
			logger.info("Waiting for mapping jobs to finish...");
			while (!clusterInterface.checkIfLoggedJobsComplete())
				Thread.sleep(5000);
	
			logger.info("OK");
	
			ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
			
			if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
				StringBuffer errMsg = new StringBuffer("The following jobs failed:\n");
				Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
				while (iteratorOverFailFiles.hasNext())
					errMsg.append(( iteratorOverFailFiles.next()).getPath() + "\n");
				logger.severe(errMsg.toString());
				throw new Exception("One or more mapping jobs failed.");
			}
			logger.info("Finished mapping jobs");
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
			throw new TaskException(e.getMessage(), e);
		}
	}


	@Override
	public boolean isDone() throws TaskException {
		checkParams();
		logger.info("Checking that mapping jobs from previous run completed successfully...");
		try {
			WholeTranscriptomeAnalyzer.checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(getMasterListOfScriptOutputFiles(), "\t", 0));
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
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
		if (remote == false) throw new IllegalArgumentException("Mapping Task only supported as Remote Task.");	
	}
	
	public void checkParams() throws TaskException {
		if (maxBytesMemoryPerJob < 100) throw new TaskException("maxBytesMemoryPerJob is too low.");
		if (memoryRequirementAdjustmentFactor < 1) throw new TaskException("memoryRequirementAdjustmentFactor cannot be less than 1");
		if (parameterTemplate == null) throw new TaskException("parameterTemplate not set");
		if (filesSplitReadsCSFasta == null) throw new TaskException("filesSplitReadsCSFasta cannot be null");
		if (filesSplitReadsCSFasta.length < 1) throw new TaskException("filesSplitReadsCSFasta cannot be empty");
		if (filesSplitReadsCSFasta.length > 2) throw new TaskException("filesSplitReadsCSFasta must have length >= 2");
		if (masksByReadSplit == null) throw new TaskException("masksByReadSplit cannot be null");
		if (masksByReadSplit.length != filesSplitReadsCSFasta.length) throw new TaskException("masksByReadSplit.length != filesSplitReadsCSFasta.length");
		if (allowedMismatches == null) throw new TaskException("allowedMismatches cannot be null");
		if (allowedMismatches.length != filesSplitReadsCSFasta.length) throw new TaskException("allowedMismatches.length != filesSplitReadsCSFasta.length");
		if (filesReferenceFastaPartitions == null) throw new TaskException("filesReferenceFastaPartitions cannot be null");
		if (filesReferenceFastaPartitions.length < 1) throw new TaskException("filesReferenceFastaPartitions.length cannot be less than 1");
		if (numberOfReads < 1) throw new TaskException("numberOfReads cannot be less than 1");
		if (applyValidAdjacentRules == null) throw new TaskException("applyValidAdjacentRules cannot be null");
		if (matchIubCodes == null) throw new TaskException("matchIubCodes cannot be null");
		if (maxMatchingLocations < 1) throw new TaskException("maxMatchingLocations cannot be less than 1");
		if (outputDir == null) throw new TaskException("outputDir cannot be null");
		if (outputDir.exists() == false) throw new TaskException("outputDir does not exist");
		if (outputDir.isDirectory() == false ) throw new TaskException("outputDir is not a directory");
		if (outputDir.canWrite() == false ) throw new TaskException("outputDir is not writeable");
		if (scratchDir == null) throw new TaskException("scratchDir cannot be null");
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


	public File[] getFilesSplitReadsCSFasta() {
		return filesSplitReadsCSFasta;
	}


	public void setFilesSplitReadsCSFasta(File[] filesSplitReadsCSFasta) throws TaskException {
		if (filesSplitReadsCSFasta.length > 2) throw new IllegalArgumentException("filesSplitReadsCSFasta cannot contain more than 2 elements.");
		this.filesSplitReadsCSFasta = filesSplitReadsCSFasta;
		this.splitLengths = new int[filesSplitReadsCSFasta.length];
		for (int i=0; i<filesSplitReadsCSFasta.length; i++) {
			try {
				this.splitLengths[i] = new FastaDatabase(this.filesSplitReadsCSFasta[i]).getMonoLength(10).intValue() - 1;

			} catch (Exception e) {
				throw new TaskException(e.getMessage(), e);
			}
		}
	}


	public ReadMask[] getMasksByReadSplit() {
		return masksByReadSplit;
	}


	public void setMasksByReadSplit(ReadMask[] masksByReadSplit) {
		this.masksByReadSplit = masksByReadSplit;
	}


	public int[] getAllowedMismatches() {
		return allowedMismatches;
	}


	public void setAllowedMismatches(int[] allowedMismatches) {
		this.allowedMismatches = allowedMismatches;
	}


	public File[] getFilesReferenceFastaPartitions() {
		return filesReferenceFastaPartitions;
	}


	public void setFilesReferenceFastaPartitions(
			File[] filesReferenceFastaPartitions) {
		this.filesReferenceFastaPartitions = filesReferenceFastaPartitions;
	}


	public int getNumberOfReads() {
		return numberOfReads;
	}


	public void setNumberOfReads(int numberOfReads) {
		this.numberOfReads = numberOfReads;
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


	public int getMaxMatchingLocations() {
		return maxMatchingLocations;
	}


	public void setMaxMatchingLocations(int maxMatchingLocations) {
		this.maxMatchingLocations = maxMatchingLocations;
	}


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


	public int getNumberOfReadsFilePartitions() {
		return numberOfReadsFilePartitions;
	}


	public void setNumberOfReadsFilePartitions(int numberOfReadsFilePartitions) {
		throw new UnsupportedOperationException("setNumberOfReadsFilePartitions() is not supported.");
	}


	public File getJunctionReference() {
		return junctionReference;
	}


	public void setJunctionReference(File junctionReference) {
		this.junctionReference = junctionReference;
	}


	public File getMappingAndExtensionFolder() {
		return new File(getTmpDir(), "mappingAndExtension");
	}
	
	public File getScriptDir() {
		return new File(getOutputDir(), "scripts");
	}
	
	public File getTmpDir() {
		return new File(getOutputDir(), "tmp");
	}
	
	public File getSplitSubDir(int leftOrRight) {
		return new File(getMappingAndExtensionFolder(), "reads_"+leftOrRight);
	}
	
	public File getPartitionSubDir(int leftOrRight, int partition) {
		return new File(getSplitSubDir(leftOrRight), "ref_"+partition);
	}
	
	public File getMaFile(int leftOrRight, int partition, int splitLength) {
		String nameOfMAFile = filesSplitReadsCSFasta[leftOrRight].getName() + "." + splitLength + "." + this.allowedMismatches[leftOrRight];
		if (this.applyValidAdjacentRules)
			nameOfMAFile += ".adj_valid";
		nameOfMAFile += ".0_" + (numberOfReads-1) + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE;
		return new File(getPartitionSubDir(leftOrRight, partition), nameOfMAFile);
	}
	
	public File[] getLeftMaFiles() {
		File[] files = new File[this.filesReferenceFastaPartitions.length];
		for (int i=0; i<this.filesReferenceFastaPartitions.length; i++) {
			files[i] = getMaFile(0, i, this.splitLengths[0]);
		}
		return files;
	}
	
	public File[] getRightMaFiles() {
		if (this.getFilesSplitReadsCSFasta().length < 2) return new File[0];
		File[] files = new File[this.filesReferenceFastaPartitions.length];
		for (int i=0; i<this.filesReferenceFastaPartitions.length; i++) {
			files[i] = getMaFile(1, i, this.splitLengths[1]);
		}
		return files;
	}
	
	public File getJunctionMaFile() {
		if (this.junctionReference == null) return null;
		String junctionMaFileName = filesSplitReadsCSFasta[0].getName()+"."+splitLengths[0]+"."+allowedMismatches[0];
		if (this.applyValidAdjacentRules)
			junctionMaFileName += ".adj_valid";
		junctionMaFileName += "."+0+"_"+(numberOfReads-1)+WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE;
		return new File(getJunctionOutputDir(), junctionMaFileName);
	}
	
	public File getJunctionOutputDir() {
		if (this.junctionReference == null) return null;
		return new File(getMappingAndExtensionFolder(), "reads_0/junctions");
	}
	
	public File getMasterJobSubmissionScript() {
		return new File (getScriptDir(), "mapping.master.script_submission.sh");
	}
	
	public File getMasterJobRemovalScript() {
		return new File(getScriptDir(), "mapping.master.script_removal.sh");
	}
	
	public File getMasterListOfScriptOutputFiles() {
		return new File(getOutputDir(), "tmp/mapping.list_of_script_output_files.out");
	}
}
