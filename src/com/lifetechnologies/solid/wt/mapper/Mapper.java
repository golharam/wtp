package com.lifetechnologies.solid.wt.mapper;

import java.io.File;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.util.FileUtils;

/**
 * main API entrypoint for mapping reads to a reference.
 * @author mullermw
 *
 */
public class Mapper {

	public Mapper() {}
	
	private Logger logger = Logger.getLogger(Mapper.class.getSimpleName());
	
	private String schedulingEnvironment;
	private String nameOfQueue;
	private long maxMemoryPerJob;
	private double memoryRequirementAdjustmentFactor;
	private String schedulerResourceRequirements;
	private String additionalSchedulerOptions;
	private MappingStep skipTo;
	private FilteringMode filteringMode;
	private int readLength;
	private ReadMask readMask;
	private int maxMappingLocations;
	private boolean validAdjacentRules;
	private boolean matchIub;
	private int leftLength;
	private int rightLength;
	private int leftMismatches;
	private int rightMismatches;
	private int filterMismatches;
	private int minScore;
	private int scoreGap;
//	private int matchScore = 1;
//	private int mismatchPenalty = -1;
	private int knownJunctionPenalty;
	private int putativeJunctionPenalty;
	private File filterReference;
	private File referenceFile;
	private File exonReference;
	private File csfastaFile;
	private File scratchDir;
	private File outputDir;
	private boolean deleteIntermediateFiles = false;
	
	public static enum MappingStep {
		splitting, filtering, ref_partitioning, mapping, extension, merge;
		
		public static MappingStep parseMappingStep(String s) {
			return MappingStep.valueOf(s.toLowerCase());
		}
	}
	
	public String getSchedulingEnvironment() {
		return schedulingEnvironment;
	}

	public void setSchedulingEnvironment(String schedulingEnvironment) {
		this.schedulingEnvironment = schedulingEnvironment;
	}

	public String getNameOfQueue() {
		return nameOfQueue;
	}

	public void setNameOfQueue(String nameOfQueue) {
		this.nameOfQueue = nameOfQueue;
	}

	public long getMaxMemoryPerJob() {
		return maxMemoryPerJob;
	}

	public void setMaxMemoryPerJob(long maxMemoryPerJob) {
		this.maxMemoryPerJob = maxMemoryPerJob;
	}

	public double getMemoryRequirementAdjustmentFactor() {
		return memoryRequirementAdjustmentFactor;
	}

	public void setMemoryRequirementAdjustmentFactor(double memoryRequirementAdjustmentFactor) {
		this.memoryRequirementAdjustmentFactor = memoryRequirementAdjustmentFactor;
	}

	public String getSchedulerResourceRequirements() {
		return schedulerResourceRequirements;
	}

	public void setSchedulerResourceRequirements(
			String schedulerResourceRequirements) {
		this.schedulerResourceRequirements = schedulerResourceRequirements;
	}

	public String getAdditionalSchedulerOptions() {
		return additionalSchedulerOptions;
	}

	public void setAdditionalSchedulerOptions(String additionalSchedulerOptions) {
		this.additionalSchedulerOptions = additionalSchedulerOptions;
	}

	public MappingStep getSkipTo() {
		return skipTo;
	}

	public void setSkipTo(MappingStep skipTo) {
		this.skipTo = skipTo;
	}

	public FilteringMode getFilteringMode() {
		return filteringMode;
	}

	public void setFilteringMode(FilteringMode filteringMode) {
		this.filteringMode = filteringMode;
	}

	public int getReadLength() {
		return readLength;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public ReadMask getReadMask() {
		return readMask;
	}

	public void setReadMask(ReadMask readMask) {
		this.readMask = readMask;
	}

	public int getMaxMappingLocations() {
		return maxMappingLocations;
	}

	public void setMaxMappingLocations(int maxMappingLocations) {
		this.maxMappingLocations = maxMappingLocations;
	}

	public boolean isValidAdjacentRules() {
		return validAdjacentRules;
	}

	public void setValidAdjacentRules(boolean validAdjacentRules) {
		this.validAdjacentRules = validAdjacentRules;
	}

	public boolean isMatchIub() {
		return matchIub;
	}

	public void setMatchIub(boolean matchIub) {
		this.matchIub = matchIub;
	}

	public int getLeftLength() {
		return leftLength;
	}

	public void setLeftLength(int leftLength) {
		this.leftLength = leftLength;
	}

	public int getRightLength() {
		return rightLength;
	}

	public void setRightLength(int rightLength) {
		this.rightLength = rightLength;
	}

	public int getLeftMismatches() {
		return leftMismatches;
	}

	public void setLeftMismatches(int leftMismatches) {
		this.leftMismatches = leftMismatches;
	}

	public int getRightMismatches() {
		return rightMismatches;
	}

	public void setRightMismatches(int rightMismatches) {
		this.rightMismatches = rightMismatches;
	}

	public int getFilterMismatches() {
		return filterMismatches;
	}

	public void setFilterMismatches(int filterMismatches) {
		this.filterMismatches = filterMismatches;
	}

	public int getMinScore() {
		return minScore;
	}

	public void setMinScore(int minScore) {
		this.minScore = minScore;
	}

	public int getScoreGap() {
		return scoreGap;
	}

	public void setScoreGap(int scoreGap) {
		this.scoreGap = scoreGap;
	}

	public File getFilterReference() {
		return filterReference;
	}

	public void setFilterReference(File filterReference) {
		this.filterReference = filterReference;
	}

	public File getReferenceFile() {
		return referenceFile;
	}

	public void setReferenceFile(File referenceFile) {
		this.referenceFile = referenceFile;
	}

	public File getExonReference() {
		return exonReference;
	}

	public void setExonReference(File exonReference) {
		this.exonReference = exonReference;
	}

	public File getCsfastaFile() {
		return csfastaFile;
	}

	public void setCsfastaFile(File csfastaFile) {
		this.csfastaFile = csfastaFile;
	}

	public File getScratchDir() {
		return scratchDir;
	}

	public void setScratchDir(File scratchDir) {
		this.scratchDir = scratchDir;
	}

	public File getOutputDir() {
		return outputDir;
	}
	
	public File getTmpDir() {
		return new File(outputDir, "tmp");
	}
	
	public File getScriptDir() {
		return new File(outputDir, "scripts");
	}
	
	public File getResultDir() {
		return new File(outputDir, "output");
	}

	public void setOutputDir(File outputDir) {
		this.outputDir = outputDir;
	}
//	
//	public int getMatchScore() {
//		return matchScore;
//	}
//
//	public void setMatchScore(int matchScore) {
//		this.matchScore = matchScore;
//	}
//
//	public int getMismatchPenalty() {
//		return mismatchPenalty;
//	}
//
//	public void setMismatchPenalty(int mismatchPenalty) {
//		this.mismatchPenalty = mismatchPenalty;
//	}

	public int getKnownJunctionPenalty() {
		return knownJunctionPenalty;
	}

	public void setKnownJunctionPenalty(int knownJunctionPenalty) {
		this.knownJunctionPenalty = knownJunctionPenalty;
	}

	public int getPutativeJunctionPenalty() {
		return putativeJunctionPenalty;
	}

	public void setPutativeJunctionPenalty(int putativeJunctionPenalty) {
		this.putativeJunctionPenalty = putativeJunctionPenalty;
	}

	public boolean isDeleteIntermediateFiles() {
		return deleteIntermediateFiles;
	}

	public void setDeleteIntermediateFiles(boolean deleteIntermediateFiles) {
		this.deleteIntermediateFiles = deleteIntermediateFiles;
	}

	public void mapReads() throws Exception {
		
		if (!outputDir.exists()) outputDir.mkdir();
		File tmpDir = getTmpDir();
		if (!tmpDir.exists()) tmpDir.mkdir();
		File scriptDir = getScriptDir();
		if (!scriptDir.exists()) scriptDir.mkdir();
		File resultDir = getResultDir();
		if (!resultDir.exists()) resultDir.mkdir();
		
		JobSubmissionParameters jobTemplate = new JobSubmissionParameters();
		jobTemplate.setEnvironment(this.getSchedulingEnvironment());
		jobTemplate.setQueueName(this.getNameOfQueue());
		jobTemplate.setResourceString(this.getSchedulerResourceRequirements());
		jobTemplate.setAdditionalOptions(this.getAdditionalSchedulerOptions());
		
		SplitReadsTask splitReadsTask = new SplitReadsTask(this.csfastaFile, tmpDir, this.leftLength, this.rightLength);
		splitReadsTask.setJobSubmissionParameters(jobTemplate.clone());
		splitReadsTask.setRemote(true);
		
		boolean mappingRightSide = splitReadsTask.getSplitReadsFiles().length > 1;
		int[] mismatches = mappingRightSide ? new int[] {this.leftMismatches, this.rightMismatches} : new int[] {this.leftMismatches};
		int[] splitLengths = mappingRightSide ? new int[] {this.leftLength, this.rightLength} : new int[] {this.leftLength};
		ReferencePartitioningTask refPartitioningTask = new ReferencePartitioningTask(this.referenceFile, this.maxMemoryPerJob, tmpDir, this.readLength, this.exonReference );
		refPartitioningTask.setJobSubmissionParameters(jobTemplate.clone());
		refPartitioningTask.setRemote(true);
		
		FilteringTask filteringTask = new FilteringTask();
		filteringTask.setFilteringMode(this.filteringMode);
		filteringTask.setFilterReferenceFile(this.filterReference);
		filteringTask.setJobTemplate(jobTemplate.clone());
		filteringTask.setMaxMismatchesInReadPartsForReadFiltering(this.getFilterMismatches());
		filteringTask.setMemoryRequirementAdjustmentFactor(this.getMemoryRequirementAdjustmentFactor());
		filteringTask.setOutputDir(this.getOutputDir());
		filteringTask.setScratchDir(this.getScratchDir());
		filteringTask.setSplitReadsFiles(splitReadsTask.getSplitReadsFiles());
		ReadMask leftMask = this.readMask.subMask(0, this.leftLength);
		ReadMask rightMask = new ReadMask("0"+this.readMask.subMask(this.readLength - this.rightLength, this.readLength ).asBitString());
		ReadMask[] masks = mappingRightSide ? new ReadMask[] {leftMask, rightMask} : new ReadMask[]{leftMask};
		filteringTask.setMasksByReadSplit(masks);
		logger.fine("leftMask="+leftMask);
		logger.fine("rightMask="+rightMask);
		filteringTask.setRemote(true);
		
		MappingTask mappingTask = new MappingTask();
		mappingTask.setAllowedMismatches(mismatches); 
		mappingTask.setApplyValidAdjacentRules(this.validAdjacentRules);
		mappingTask.setMasksByReadSplit(masks);
		mappingTask.setMatchIubCodes(this.matchIub);
		mappingTask.setMaxBytesMemoryPerJob(this.maxMemoryPerJob);
		mappingTask.setMaxMatchingLocations(this.maxMappingLocations + 1);
		mappingTask.setMemoryRequirementAdjustmentFactor(this.memoryRequirementAdjustmentFactor);
		mappingTask.setOutputDir(this.getOutputDir());
		mappingTask.setParameterTemplate(jobTemplate);
		mappingTask.setScratchDir(this.scratchDir);
		if (this.exonReference != null)
			mappingTask.setJunctionReference(refPartitioningTask.getJunctionFastaFile());
		mappingTask.setRemote(true);
		
		ExtensionTask extensionTask = new ExtensionTask();
		extensionTask.setApplyValidAdjacentRules(this.validAdjacentRules);
		extensionTask.setCompressIntermediateFiles(false);
		extensionTask.setDeleteIntermediateFiles(false);
		extensionTask.setFileFullLengthReadsCsfasta(this.csfastaFile);
		if (this.exonReference != null)
			extensionTask.setJunctionReference(refPartitioningTask.getJunctionFastaFile());
		extensionTask.setMatchIubCodes(this.matchIub);
		extensionTask.setMaxBytesMemoryPerJob(this.maxMemoryPerJob);
		extensionTask.setMemoryRequirementAdjustmentFactor(this.memoryRequirementAdjustmentFactor);
		extensionTask.setOutputDir(this.outputDir);
		extensionTask.setParameterTemplate(jobTemplate);
		extensionTask.setPositionOfLeftRead(1);
		extensionTask.setPositionOfRightRead(this.readLength - this.rightLength );
		extensionTask.setReadMask(this.readMask);
		extensionTask.setScratchDir(this.scratchDir);
		extensionTask.setRemote(true);
		
		MergeTask mergeTask = new MergeTask();
		mergeTask.setCompressIntermediateFiles(false);
		mergeTask.setCountsMaxErrorsAllowedWhenMatching(mismatches);
		mergeTask.setDeleteIntermediateFiles(false);
		mergeTask.setExtensionMatchScore(extensionTask.getMatchScore());
		mergeTask.setExtensionMismatchPenalty(extensionTask.getMismatchPenalty());
		mergeTask.setFilesSplitReadsCSFasta(splitReadsTask.getSplitReadsFiles());
		mergeTask.setFilteringMode(this.filteringMode);
		mergeTask.setFullLengthReadsCsfasta(this.csfastaFile);
		mergeTask.setLengthOfReads(this.readLength);
		mergeTask.setLengthsOfSplitReads(splitLengths );
		mergeTask.setMaxMappingLocationsBeforeNotReportingRead(this.maxMappingLocations);
		mergeTask.setMaxMemorytPerJobInBytes(this.maxMemoryPerJob);
		mergeTask.setMemoryRequirementAdjustmentFactor(this.memoryRequirementAdjustmentFactor);
		mergeTask.setMinAlignmentScoreForReportingAlignment(this.minScore);
		mergeTask.setMinMappingLocationsRequiredBeforeReportingRead(1);
		mergeTask.setNumberOfReadsPerMergeJobOverride(null);//TODO expose this parameter?
		mergeTask.setOutputDir(this.outputDir);
		mergeTask.setParameterTemplate(jobTemplate);
		mergeTask.setReadMasks(masks);
		mergeTask.setScratchDir(this.scratchDir);
		mergeTask.setUniquenessGap(this.scoreGap);
		mergeTask.setKnownJunctionPenalty(this.knownJunctionPenalty);
		mergeTask.setPutativeJunctionPenalty(this.putativeJunctionPenalty);
		mergeTask.setRemote(true);
			
		FastaDatabase reads = new FastaDatabase(this.csfastaFile);
		
		int readCount = -1;
		switch (skipTo) {
		
			case splitting:
				FileUtils.deleteDir(tmpDir);
				tmpDir.mkdir();
				splitReadsTask.doTask();
				readCount = splitReadsTask.getReadCount();
				
			case ref_partitioning:
				if (!splitReadsTask.isDone()) throw new TaskException("Read splitting has not been performed with these parameters.");
				refPartitioningTask.doTask();
			
			case filtering:
				if (!splitReadsTask.isDone()) throw new TaskException("Read splitting has not been performed with these parameters.");
				if (!refPartitioningTask.isDone()) throw new TaskException("Reference partitioning has not been performed with these parameters.");
				if (readCount < 0) {
					logger.info("counting number of reads");
					readCount = reads.getNumberOfSequencesInDatabase();
				}
				filteringTask.setNumberOfReads(readCount);
				filteringTask.doTask();
			case mapping:
				if (!splitReadsTask.isDone()) throw new TaskException("Read splitting has not been performed with these parameters.");
				if (!refPartitioningTask.isDone()) throw new TaskException("Reference partitioning has not been performed with these parameters.");
				if (readCount < 0) {
					logger.info("counting number of reads");
					readCount = reads.getNumberOfSequencesInDatabase();
				}
				filteringTask.setNumberOfReads(readCount);
				if (!filteringTask.isDone()) throw new TaskException("Filtering has not been performed with these parameters.");
				mappingTask.setNumberOfReads(readCount);
				mappingTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				mappingTask.setFilesSplitReadsCSFasta(filteringTask.getFilteredSplitReadsFiles());
				mappingTask.doTask();
				
			case extension:
				if (!splitReadsTask.isDone()) throw new TaskException("Read splitting has not been performed with these parameters.");
				if (!refPartitioningTask.isDone()) throw new TaskException("Reference partitioning has not been performed with these parameters.");
				if (readCount < 0) {
					logger.info("Counting number of reads");
					readCount = reads.getNumberOfSequencesInDatabase();
				}
				filteringTask.setNumberOfReads(readCount);
				if (!filteringTask.isDone()) throw new TaskException("Filtering has not been performed with these parameters.");
				mappingTask.setNumberOfReads(readCount);
				mappingTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				mappingTask.setFilesSplitReadsCSFasta(filteringTask.getFilteredSplitReadsFiles());
				if (!mappingTask.isDone()) throw new TaskException("Mapping has not been performed with these parameters.");
				
				extensionTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				extensionTask.setMaFileJunctions(mappingTask.getJunctionMaFile());
				extensionTask.setMaFilesLeft(mappingTask.getLeftMaFiles());
				extensionTask.setMaFilesRight(mappingTask.getRightMaFiles());
				extensionTask.doTask();
			case merge:
				if (!splitReadsTask.isDone()) throw new TaskException("Read splitting has not been performed with these parameters.");
				if (!refPartitioningTask.isDone()) throw new TaskException("Reference partitioning has not been performed with these parameters.");
				if (readCount < 0) {
					logger.info("Counting number of reads");
					readCount = reads.getNumberOfSequencesInDatabase();
				}
				filteringTask.setNumberOfReads(readCount);
				if (!filteringTask.isDone()) throw new TaskException("Filtering has not been performed with these parameters.");
				mappingTask.setNumberOfReads(readCount);
				mappingTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				mappingTask.setFilesSplitReadsCSFasta(filteringTask.getFilteredSplitReadsFiles());
				if (!mappingTask.isDone()) throw new TaskException("Mapping has not been performed with these parameters.");
				extensionTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				extensionTask.setMaFileJunctions(mappingTask.getJunctionMaFile());
				extensionTask.setMaFilesLeft(mappingTask.getLeftMaFiles());
				extensionTask.setMaFilesRight(mappingTask.getRightMaFiles());
				if (!extensionTask.isDone()) throw new TaskException("Extension has not been performed with these parameters.");
				mergeTask.setFilesReferenceFastaPartitions(refPartitioningTask.getPartitions());
				mergeTask.setListsOfFullRefSequenceNumbersOrderedForEachReferencePartition(refPartitioningTask.getListsOfFullRefSequenceNumbersOrderedForEachReferencePartition());
				mergeTask.setNumberOfReads(readCount);
				mergeTask.setJunctionMaxFile(extensionTask.getJunctionMaxFile());
				mergeTask.doTask();
				if (this.deleteIntermediateFiles) FileUtils.deleteDir(getTmpDir());
			}
	}
}
