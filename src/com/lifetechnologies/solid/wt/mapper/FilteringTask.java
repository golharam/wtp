package com.lifetechnologies.solid.wt.mapper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.CompressionUtilities;
import com.lifetechnologies.solid.wt.FastaDatabase;
import com.lifetechnologies.solid.wt.ReadMapper;
import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.TextFileUtilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.mapper.FilteringMode;
import com.lifetechnologies.util.FileHeader;
import com.lifetechnologies.util.FileUtils;

/**
 * The task of identifying reads that map to a filter reference.
 * @author mullermw
 *
 */
public class FilteringTask implements Task {
	
	private static Logger logger = Logger.getLogger(FilteringTask.class.toString());
	
	/* inputs */
	private File filterReferenceFile;
	private File outputDir;
	private JobSubmissionParameters jobTemplate;
	private double memoryRequirementAdjustmentFactor = 1;
	private File[] splitReadsFiles;
	private int maxMismatches;
	private ReadMask[] masksByReadSplit;
	private int numberOfReads;
	final private Boolean applyValidAdjacentRules = false;
	private File scratchDir;
	private FilteringMode filteringMode;
	
	private Boolean compressIntermediateFiles = false;
	private Boolean deleteIntermediateFiles = false;
	
	private int numberOfReadsFilePartitions = 1;
	
	/* outputs */
	private int[] splitLengths = new int[2];
	@Override
	public void doTask() throws TaskException {
		try {
			checkParams();
			logger.info("Started read filtering at " + new Date(System.currentTimeMillis()));
	        logger.info("Reads will be filtered against the following fasta file: " + filterReferenceFile.getPath());
	
	        FileUtils.deleteDir(getFilteringDir());
	        getFilterReportFile().delete();
	        getMasterFilteringJobRemovalScript().delete();
	        getMasterFilteringJobSubmissionScript().delete();
	        getMasterFilteringListOfScriptOutputFiles().delete();
	        for (File file : getTmpDir().listFiles(new FilenameFilter() {
	        	@Override
	        	public boolean accept(File dir, String name) {
	        		return name.endsWith(WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	        	}
	        })) {
	        	file.delete();
	        }
	        
	        ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();
	        
	        long memoryRequiredForMapReads = ReadMapper.calculateMemoryRequiredByMapReadsInBytes(new FastaDatabase(filterReferenceFile).getTotalLengthOfSequence());
	        memoryRequiredForMapReads = (long)Math.ceil(memoryRequiredForMapReads * memoryRequirementAdjustmentFactor);
	        
	        JobSubmissionParameters params = jobTemplate.clone();
	        params.setMemoryRequirement(memoryRequiredForMapReads);
	        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
	
	        File[] filesMAByReadSequenceSplit = new File[splitReadsFiles.length];
	        File[] filesFilteredReadsByReadSequenceSplit = new File[splitReadsFiles.length];
	
	        for (int indexSplitReadsFile = 0; indexSplitReadsFile < splitReadsFiles.length; indexSplitReadsFile++) {
	
	            File folderForFilterIO = getSplitFolder(indexSplitReadsFile);
	            if (!folderForFilterIO.exists())
	                folderForFilterIO.mkdirs();
	
	            Map<String, FileHeader.Value> header = FileHeader.parseHeader(splitReadsFiles[indexSplitReadsFile]);
	            int splitLength = header.get(FileHeader.SPLIT_LENGTH_KEY).intValue();
	            if (header.get(FileHeader.SPLIT_KEY).toString().trim().toLowerCase().equalsIgnoreCase("right")) splitLength += 1;
	            this.splitLengths[indexSplitReadsFile] = splitLength;
	            ReadMapper mapper = new ReadMapper(WholeTranscriptomeAnalyzer.getMapreadsExe(), filterReferenceFile, splitReadsFiles[indexSplitReadsFile],
	                                               splitLength, masksByReadSplit[indexSplitReadsFile].asBitString());
	
	
	            int numberOfReadsPerJob = (int)Math.ceil((double)numberOfReads / numberOfReadsFilePartitions);
	            for (int indexOfRead = 0; indexOfRead < numberOfReads; indexOfRead += numberOfReadsPerJob) {
	
	                int maxReadIndexInCurrentRange = Math.min(indexOfRead + numberOfReadsPerJob, numberOfReads) -1;
	
	                filesMAByReadSequenceSplit[indexSplitReadsFile] = getMaFile(indexSplitReadsFile, indexOfRead, maxReadIndexInCurrentRange);
	
	                mapper.startMapReadsJob(true, indexOfRead, maxReadIndexInCurrentRange,
	                                        maxMismatches, this.applyValidAdjacentRules, true, 1,
	                                        clusterInterface, "wt_filter.",
	                                        scratchDir, WholeTranscriptomeAnalyzer.getSchemaDir(), folderForFilterIO, filesMAByReadSequenceSplit[indexSplitReadsFile], true);
	                 
	                listOfFilesThatAreNoLongerNeeded.add(filesMAByReadSequenceSplit[indexSplitReadsFile]);
	            }
	
	            filesFilteredReadsByReadSequenceSplit[indexSplitReadsFile] = getFilteredSplitReadsFile(indexSplitReadsFile);
	        }
	
	        File masterSubmissionScript = getMasterFilteringJobSubmissionScript();
	        clusterInterface.writeMasterJobSubmissionFileFromLog(masterSubmissionScript);
	        File masterRemovalScriptFile = getMasterFilteringJobRemovalScript();
	        clusterInterface.writeMasterJobRemovalFileFromLog(masterRemovalScriptFile);
	        File masterListOfScriptOutputFiles = getMasterFilteringListOfScriptOutputFiles();
	        clusterInterface.writeMasterListOfJobOutputFilesFromLog(masterListOfScriptOutputFiles);
	
	        masterSubmissionScript.setExecutable(true, true);
	        masterRemovalScriptFile.setExecutable(true, true);
	
	        // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
	        logger.info("Waiting for filtering jobs to finish...");
	        while (!clusterInterface.checkIfLoggedJobsComplete())
	            Thread.sleep(5000);
	
	        logger.info("OK. Detected jobs complete.");
	
	        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
	        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
	        	StringBuffer msg = new StringBuffer("The following jobs failed:\n");
	        	for (File file : listOfScriptOutputFilesForFailedJobs)
	                msg.append(file.getPath()+"\n");
	        	logger.severe(msg.toString());
	            throw new Exception("One or more filtering jobs failed.");
	        }
	
	        logger.info("Tagging filtered reads and generating a filter report.");
	
	        processFilteredReadsAndWriteReport(true,
	                                            maxMismatches,
	                                            null,
	                                            filterReferenceFile,
	                                            filesMAByReadSequenceSplit,
	                                            filesFilteredReadsByReadSequenceSplit,
	                                            getFilterReportFile());
	
	        //if (filteringMode != FilteringMode.OFF)
	            //for (int indexSplitReadsFile = 0; indexSplitReadsFile < splitReadsFiles.length; indexSplitReadsFile++)
	            //    splitReadsFiles[indexSplitReadsFile] = filesFilteredReadsByReadSequenceSplit[indexSplitReadsFile];
	
	
	        if (deleteIntermediateFiles || compressIntermediateFiles) {
	            if (compressIntermediateFiles && deleteIntermediateFiles)  
	            	logger.info("Compressing and deleting intermediate files...");
	            else if (compressIntermediateFiles)   logger.info("Compressing intermediate files...");
	            else    logger.info("Deleting intermediate files...");
	            Iterator<File> iteratorFilesNoLongerNeeded = listOfFilesThatAreNoLongerNeeded.iterator();
	            while (iteratorFilesNoLongerNeeded.hasNext()) {
	                File file = iteratorFilesNoLongerNeeded.next();
	                if (compressIntermediateFiles)
	                    CompressionUtilities.compressFile(file, new File(file.getPath() + WholeTranscriptomeAnalyzer.EXTENSION_FOR_COMPRESSED_FILE));
	                if (deleteIntermediateFiles)
	                    file.delete();
	            }
	            logger.info("OK. Finished deleting or compressing intermediate files.");
	        }
		} catch (Exception e) {
			logger.log(Level.SEVERE, "Failure in Filtering", e);
			throw new TaskException("Failure in Filtering", e);
		}
		logger.info("Finished read filtering");
    }
	
    
    /**
    *
    * @param justMarkAsFilteredAndDontActuallyFilter if true, reads with too few/many mappings will still be reported, but their headings will be marked to indicate they should be filtered and mappings will not be reported
    * @param maxErrorsAllowedWhenMatching
    * //@param listsOfFullRefSequenceNumbersOrderedForEachReferencePartition  this data structure specifies the following mapping:  reference_partition_index -> sequence_number_in_ref_partition -> sequence_number_in_full_reference
    *   //                                                                    if not specified, then sequence numbers will be the same in the merged file as in the partitioned MA files (not recommended unless # of ref file paritions == 1).
    * @param fileOfFullLengthReadsCSFasta if this file is specified, then the sequences contained therein will be printed in the merged file.
    *                                      this file may contain more read sequences than are found in the split MA files, but must be ordered in the same fashion.

    * @throws Exception
    */
    public static void processFilteredReadsAndWriteReport(boolean justMarkAsFilteredAndDontActuallyFilter,
                                                          int maxErrorsAllowedWhenMatching,
                                                          //ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
                                                          File fileOfFullLengthReadsCSFasta,
                                                          File fileReferenceFasta,
                                                          File[] filesMAByReadSequenceSplit,
                                                          File[] filesFilteredReadsByReadSequenceSplit,
                                                          File fileReport) throws Exception {

       FastaDatabase fastaDatabaseReference = new FastaDatabase(fileReferenceFasta);
       List<String> listOfReferenceHeaders = fastaDatabaseReference.getListOfHeaders(false, false);
       // this data structure will track matches to each reference sequence stratified by read sequence split and mismatch count
       long countOfMatchesByRefSequenceByReadSplitByMismatchCount[][][] = new long[listOfReferenceHeaders.size()][filesMAByReadSequenceSplit.length][maxErrorsAllowedWhenMatching +1];
       long countOfMatchesByRefSequenceWhereAllSplitsAgree[] = new long[listOfReferenceHeaders.size()];

       BufferedReader readerFullLengthReadsCSFasta = null;
       if (fileOfFullLengthReadsCSFasta != null)
           readerFullLengthReadsCSFasta = new BufferedReader(new FileReader(fileOfFullLengthReadsCSFasta));

       BufferedReader readersMAByReadSplit[] = new BufferedReader[filesMAByReadSequenceSplit.length];
       BufferedWriter writersMergedOutputByReadSplit[] = new BufferedWriter[filesFilteredReadsByReadSequenceSplit.length];
       for (int i = 0; i < writersMergedOutputByReadSplit.length; i++) {
           readersMAByReadSplit[i] = new BufferedReader(new FileReader(filesMAByReadSequenceSplit[i]));
           writersMergedOutputByReadSplit[i] = new BufferedWriter(new FileWriter(filesFilteredReadsByReadSequenceSplit[i]));
       }

       long countOfReadsProcessed = 0;
       long countOfReadsWithMappingsBySplit[] = new long[filesMAByReadSequenceSplit.length];

       // drop the comments at the top of each match file
       String lines[] = new String[readersMAByReadSplit.length];
       for (int i = 0; i < readersMAByReadSplit.length; i++) {
           lines[i] = readersMAByReadSplit[i].readLine();
           while (lines[i] != null && !lines[i].startsWith(">"))
               lines[i] = readersMAByReadSplit[i].readLine();
       }

       while (lines[0] != null) {  // each loop iteration deals with one read

           countOfReadsProcessed++;

           int indexOfFirstComma = lines[0].indexOf(',');
           StringBuffer headerWithoutMappingInfo = null;
           if (indexOfFirstComma > 0)
               headerWithoutMappingInfo = new StringBuffer(lines[0].substring(0, indexOfFirstComma));
           else
               headerWithoutMappingInfo = new StringBuffer(lines[0]);


           String sequenceOfReadFull = null;
           if (readerFullLengthReadsCSFasta != null) {
               // will replace split read sequence with the original, full read sequence
               // assumption is that the reads are ordered the same way in this file and all MA files
               // however, this file can contain a superset of reads --> therefore must check headers and skip ahead if necessary
               boolean advancedReaderToCurrentRead = false;
               while (!advancedReaderToCurrentRead) {
                   String headerOfCurrentRead = readerFullLengthReadsCSFasta.readLine();
                   if (headerOfCurrentRead.equalsIgnoreCase(headerWithoutMappingInfo.toString()))
                       advancedReaderToCurrentRead = true;
                   sequenceOfReadFull = readerFullLengthReadsCSFasta.readLine();
               }
           }

           boolean foundMatchByReadSplitLocal[] = new boolean[filesFilteredReadsByReadSequenceSplit.length];
           HashMap<Integer, boolean[][]> mapSequenceNumberToFoundMatchByReadSplitByMismatchCount = new HashMap<Integer, boolean[][]>();
           for (int indexReadSplit = 0; indexReadSplit < readersMAByReadSplit.length; indexReadSplit++) {
               String coordinatesOfMatches[] = lines[indexReadSplit].split(",");
               for (int indexOfMatch = 1; indexOfMatch < coordinatesOfMatches.length; indexOfMatch++) {
                   String tokens[] = coordinatesOfMatches[indexOfMatch].split("_");
                   int indexOfSequence = Integer.parseInt(tokens[0]) -1;
                   int numberOfMismatches = Integer.parseInt(tokens[1].split("\\.")[1]);
                   if (numberOfMismatches <= maxErrorsAllowedWhenMatching) {
                       if (mapSequenceNumberToFoundMatchByReadSplitByMismatchCount.containsKey(indexOfSequence))
                           mapSequenceNumberToFoundMatchByReadSplitByMismatchCount.get(indexOfSequence)[indexReadSplit][numberOfMismatches] = true;
                       else {
                           boolean foundMatchByReadSplitByMismatchCount[][] = new boolean[filesMAByReadSequenceSplit.length][maxErrorsAllowedWhenMatching +1];
                           foundMatchByReadSplitByMismatchCount[indexReadSplit][numberOfMismatches] = true;
                           mapSequenceNumberToFoundMatchByReadSplitByMismatchCount.put(indexOfSequence, foundMatchByReadSplitByMismatchCount);
                       }
                       foundMatchByReadSplitLocal[indexReadSplit] = true;

                   }
               }
               String sequenceOfReadSplit = readersMAByReadSplit[indexReadSplit].readLine();
               lines[indexReadSplit] = readersMAByReadSplit[indexReadSplit].readLine();

               if (!foundMatchByReadSplitLocal[indexReadSplit] || justMarkAsFilteredAndDontActuallyFilter) {
                   writersMergedOutputByReadSplit[indexReadSplit].write(headerWithoutMappingInfo.toString());
                   if (foundMatchByReadSplitLocal[indexReadSplit])
                       writersMergedOutputByReadSplit[indexReadSplit].write(WholeTranscriptomeAnalyzer.SUFFIX_FOR_MARKING_READS_AS_FILTERED);
                   writersMergedOutputByReadSplit[indexReadSplit].newLine();
                   if (sequenceOfReadFull == null)
                       writersMergedOutputByReadSplit[indexReadSplit].write(sequenceOfReadSplit);
                   else
                       writersMergedOutputByReadSplit[indexReadSplit].write(sequenceOfReadFull);
                   writersMergedOutputByReadSplit[indexReadSplit].newLine();
               }

           }

           Iterator<Integer> iteratorOverSequencesThatWereMatched = mapSequenceNumberToFoundMatchByReadSplitByMismatchCount.keySet().iterator();
           while (iteratorOverSequencesThatWereMatched.hasNext()) {
               int indexOfSequence = iteratorOverSequencesThatWereMatched.next();
               boolean foundMatchByReadSplitByMismatchCount[][] = mapSequenceNumberToFoundMatchByReadSplitByMismatchCount.get(indexOfSequence);
               boolean foundMappingForAllSplits = true;
               for (int indexReadSplit = 0; indexReadSplit < foundMatchByReadSplitByMismatchCount.length; indexReadSplit++) {
                   boolean foundMappingForThisSplit = false;
                   for (int indexMismatchCount = 0; indexMismatchCount < foundMatchByReadSplitByMismatchCount[indexReadSplit].length; indexMismatchCount++)
                       if (foundMatchByReadSplitByMismatchCount[indexReadSplit][indexMismatchCount]) {
                           countOfMatchesByRefSequenceByReadSplitByMismatchCount[indexOfSequence][indexReadSplit][indexMismatchCount]++;
                           foundMappingForThisSplit = true;
                       }
                   foundMappingForAllSplits &= foundMappingForThisSplit;
               }
               if (foundMappingForAllSplits)
                   countOfMatchesByRefSequenceWhereAllSplitsAgree[indexOfSequence]++;
           }

           for (int indexReadSplit = 0; indexReadSplit < readersMAByReadSplit.length; indexReadSplit++) {
               if (foundMatchByReadSplitLocal[indexReadSplit])
                   countOfReadsWithMappingsBySplit[indexReadSplit]++;
           }


       }


       if (readerFullLengthReadsCSFasta != null)
           readerFullLengthReadsCSFasta.close();
       for (int i = 0; i < readersMAByReadSplit.length; i++) {
           readersMAByReadSplit[i].close();
           writersMergedOutputByReadSplit[i].close();
       }

       BufferedWriter writerReport = new BufferedWriter(new FileWriter(fileReport));

       logger.info("Writing filter report " + fileReport.getPath());

       DecimalFormat formatter = new DecimalFormat("0");

       logger.info("countOfReadsProcessed:\t" + formatter.format(countOfReadsProcessed));
       writerReport.write("countOfReadsProcessed:\t" + formatter.format(countOfReadsProcessed));
       writerReport.newLine();
       for (int i = 0; i < countOfReadsWithMappingsBySplit.length; i++) {
           String percentage = " (" + (new DecimalFormat("##0.0")).format(100.0 * countOfReadsWithMappingsBySplit[i] / countOfReadsProcessed) + "%)";
           logger.info("countOfReadsWithMappingsBySplit_" + i + ":\t" + formatter.format(countOfReadsWithMappingsBySplit[i]) + percentage);
           writerReport.write("countOfReadsWithMappingsBySplit_" + i + ":\t" + formatter.format(countOfReadsWithMappingsBySplit[i]) + percentage);
           writerReport.newLine();
       }
       writerReport.newLine();

       for (int indexReadSplit = 0; indexReadSplit < filesMAByReadSequenceSplit.length; indexReadSplit++)
           for (int indexMismatchCount = 0; indexMismatchCount < maxErrorsAllowedWhenMatching +2; indexMismatchCount++)
               writerReport.write("\treadSplit_" + indexReadSplit);
       writerReport.write("\treadSplitsThatAgree");

       writerReport.newLine();
       writerReport.write("sequenceHeader");
       for (int indexReadSplit = 0; indexReadSplit < filesMAByReadSequenceSplit.length; indexReadSplit++) {
           for (int indexMismatchCount = 0; indexMismatchCount < maxErrorsAllowedWhenMatching +1; indexMismatchCount++)
               writerReport.write("\t" + indexMismatchCount + "_mismatches");
           writerReport.write("\ttotal");
       }
       writerReport.write("\ttotal");
       writerReport.newLine();

       for (int indexRefSequence = 0; indexRefSequence < countOfMatchesByRefSequenceByReadSplitByMismatchCount.length; indexRefSequence++) {

           writerReport.write(listOfReferenceHeaders.get(indexRefSequence).replaceAll("\t", "_"));
           for (int indexReadSplit = 0; indexReadSplit < countOfMatchesByRefSequenceByReadSplitByMismatchCount[indexRefSequence].length; indexReadSplit++) {
               long sumOfCountsForReadSplit = 0;
               for (int indexMismatchCount = 0; indexMismatchCount < countOfMatchesByRefSequenceByReadSplitByMismatchCount[indexRefSequence][indexReadSplit].length; indexMismatchCount++) {
                   writerReport.write("\t" + formatter.format(countOfMatchesByRefSequenceByReadSplitByMismatchCount[indexRefSequence][indexReadSplit][indexMismatchCount]));
                   sumOfCountsForReadSplit += countOfMatchesByRefSequenceByReadSplitByMismatchCount[indexRefSequence][indexReadSplit][indexMismatchCount];
               }
               writerReport.write("\t" + formatter.format(sumOfCountsForReadSplit));
           }
           writerReport.write("\t" + formatter.format(countOfMatchesByRefSequenceWhereAllSplitsAgree[indexRefSequence]));
           writerReport.newLine();
       }

       writerReport.close();
   }
    
	@Override
	public boolean isDone() throws TaskException {
		checkParams();
        logger.info("Checking that filtering jobs from previous run completed successfully...");
        if (this.filteringMode == FilteringMode.OFF) return true;
        try {
        	WholeTranscriptomeAnalyzer.checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(getMasterFilteringListOfScriptOutputFiles(), "\t", 0));
        } catch (Exception e) {
        	logger.log(Level.SEVERE, e.getMessage(), e);
        	return false;
        }
        logger.info(" OK");
        File[] filesSplitReadsCSFasta = getSplitReadsFiles();
        for (int indexReadsFile = 0; indexReadsFile < filesSplitReadsCSFasta.length; indexReadsFile++) {
            File fileMergedFilteredOutput = new File(filesSplitReadsCSFasta[indexReadsFile].getPath().replaceAll(WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE + "$", WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE));
            logger.info("Checking that post-filter reads file exist for " + filesSplitReadsCSFasta[indexReadsFile]+"->"+fileMergedFilteredOutput);
            if (!fileMergedFilteredOutput.exists()) {
                logger.severe("Could not find filtered reads file: " + fileMergedFilteredOutput + ".  Please re-run with RUN_READ_FILTERING on.");
            	return false;
            }
        }
        return true;
	}
	
	@Override
	public boolean isRemote() {
		return true;
	}
	
	@Override
	public void setRemote(boolean remote) {
		if (remote == false) throw new IllegalArgumentException("FilteringTask is only supported as a remote task.");
	}
	
	private void checkParams() throws TaskException {
		if (filterReferenceFile == null) throw new TaskException("filterReferenceFile cannot be null.");
		if (outputDir == null) throw new TaskException("outputDir cannot be null.");
		if (jobTemplate == null) throw new TaskException("jobTemplate cannot be null.");
		if (memoryRequirementAdjustmentFactor < 1) throw new TaskException("memoryRequirementAdjustmentFactor must be greater than 1.");
		if (splitReadsFiles == null ) throw new TaskException("splitReadsFiles cannot be null");
		if (splitReadsFiles.length < 1 || splitReadsFiles.length > 2) throw new TaskException("splitReadsFiles must contain one or two files.");
		if (masksByReadSplit.length < splitReadsFiles.length) throw new TaskException("masksByReadSplit doesn't match splitReadsFiles");
		if (numberOfReads < 1) throw new TaskException("numberOfReads cannot be less than 1.");
		if (applyValidAdjacentRules == null) throw new TaskException("applyValidAdjacentRules not set.");
		if (maxMismatches < 0) throw new TaskException("maxMismatches cannot be less than 1.");
		if (scratchDir == null) throw new TaskException("scratchDif cannot be null.");
		if (filteringMode == null) throw new TaskException("filteringMode cannot be null.");
	}
	
	public File getFilterReferenceFile() {
		return filterReferenceFile;
	}


	public void setFilterReferenceFile(File filterReferenceFile) {
		this.filterReferenceFile = filterReferenceFile;
	}


	public JobSubmissionParameters getJobTemplate() {
		return jobTemplate;
	}


	public void setJobTemplate(JobSubmissionParameters jobTemplate) {
		this.jobTemplate = jobTemplate;
	}


	public double getMemoryRequirementAdjustmentFactor() {
		return memoryRequirementAdjustmentFactor;
	}


	public void setMemoryRequirementAdjustmentFactor(
			double memoryRequirementAdjustmentFactor) {
		this.memoryRequirementAdjustmentFactor = memoryRequirementAdjustmentFactor;
	}


	public File[] getSplitReadsFiles() {
		return splitReadsFiles;
	}


	public void setSplitReadsFiles(File[] splitReadsFiles) {
		this.splitReadsFiles = splitReadsFiles;
	}


	public ReadMask[] getMasksByReadSplit() {
		return masksByReadSplit;
	}


	public void setMasksByReadSplit(ReadMask[] masksByReadSplit) {
		this.masksByReadSplit = masksByReadSplit;
	}

	public int getNumberOfReads() {
		return numberOfReads;
	}


	public void setNumberOfReads(int numberOfReads) {
		this.numberOfReads = numberOfReads;
	}


	public Boolean getValidAdjacentMismatchesCountAsOneForMappingAndExtension() {
		return applyValidAdjacentRules;
	}


	public void setApplyValidAdjacentRules(
			Boolean applyValidAdjacentRules) {
		throw new UnsupportedOperationException("operation not supported");
		//this.applyValidAdjacentRules = applyValidAdjacentRules;
	}


	public int getMaxMismatchesInReadPartsForReadFiltering() {
		return maxMismatches;
	}


	public void setMaxMismatchesInReadPartsForReadFiltering(
			int maxMismatchesInReadPartsForReadFiltering) {
		this.maxMismatches = maxMismatchesInReadPartsForReadFiltering;
	}


	public File getScratchDir() {
		return scratchDir;
	}


	public void setScratchDir(File scratchDir) {
		this.scratchDir = scratchDir;
	}


	public FilteringMode getFilteringMode() {
		return filteringMode;
	}


	public void setFilteringMode(FilteringMode filteringMode) {
		this.filteringMode = filteringMode;
	}


	public Boolean getCompressIntermediateFiles() {
		return compressIntermediateFiles;
	}


	public void setCompressIntermediateFiles(Boolean compressIntermediateFiles) {
		this.compressIntermediateFiles = compressIntermediateFiles;
	}


	public Boolean getDeleteIntermediateFiles() {
		return deleteIntermediateFiles;
	}


	public void setDeleteIntermediateFiles(Boolean deleteIntermediateFiles) {
		this.deleteIntermediateFiles = deleteIntermediateFiles;
	}


	public int getNumberOfReadsFilePartitions() {
		return numberOfReadsFilePartitions;
	}


	public void setNumberOfReadsFilePartitions(int numberOfReadsFilePartitions) {
		throw new UnsupportedOperationException("setNumberOfReadsFilePartitions() is not supported");
	}


	public void setOutputDir(File outputDir) {
		this.outputDir = outputDir;
	}


	public File getOutputDir() {
		return outputDir;
	}
	
	public File getScriptDir() {
		return new File(getOutputDir(), "scripts");
	}
	
	public File getTmpDir() {
		return new File(getOutputDir(), "tmp");
	}
	
	public File getFilteringDir() {
		return new File (getOutputDir(), "tmp/filtering");
	}
	
	public File getMasterFilteringJobSubmissionScript() {
		return new File (getScriptDir(), "filtering.master.script_submission.sh");
	}
	
	public File getMasterFilteringJobRemovalScript() {
		return new File(getScriptDir(), "filtering.master.script_removal.sh");
	}
	
	public File getMasterFilteringListOfScriptOutputFiles() {
		return new File(getOutputDir(), "tmp/filtering.list_of_script_output_files.out");
	}
	
	public File getFilterReportFile() {
		return new File(getOutputDir(), "output/filterReport.txt");
	}
	
	public File getFilteredSplitReadsFile(int leftOrRight) {
		//if (this.filteringMode == FilteringMode.OFF) return splitReadsFiles[leftOrRight];
		return new File(getTmpDir(), splitReadsFiles[leftOrRight].getName().replaceAll(WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE + "$", WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE));
	}
	
	public File[] getFilteredSplitReadsFiles() {
		File[] arr = new File[this.splitReadsFiles.length];
		for (int i=0; i<arr.length; i++) {
			arr[i] = getFilteredSplitReadsFile(i);
		}
		return arr;
	}
	
	public File getSplitFolder(int leftOrRight) {
		return new File(getFilteringDir(), "reads_" + leftOrRight + "/ref_0/");
	}
	
	public File getMaFile(int leftOrRight, int indexOfRead, int maxReadIndexInCurrentRange) {
        String nameOfMAFile = this.splitReadsFiles[leftOrRight].getName() + "." + this.splitLengths[leftOrRight] + "." + maxMismatches;
        if (this.applyValidAdjacentRules) nameOfMAFile += ".adj_valid";
        	nameOfMAFile += "." + indexOfRead + "_" + maxReadIndexInCurrentRange + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_OUTPUT_FILE;
        return new File(getSplitFolder(leftOrRight), nameOfMAFile);
	}
}
