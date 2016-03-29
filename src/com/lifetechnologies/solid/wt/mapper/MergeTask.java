package com.lifetechnologies.solid.wt.mapper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.CompressionUtilities;
import com.lifetechnologies.solid.wt.ExtendedReadMapping;
import com.lifetechnologies.solid.wt.ExtendedReadMappingSetForRead;
import com.lifetechnologies.solid.wt.ExtendedReadMappingsMerger;
import com.lifetechnologies.solid.wt.ExtendedReadMappingsReport;
import com.lifetechnologies.solid.wt.ReadMappingComparator;
import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.config.WTConfigUtils;
import com.lifetechnologies.util.FileUtils;
import com.lifetechnologies.util.LineUtils;

/**
 * The task of merging the results of many extension tasks into a single result.
 * Prepares reports.
 * Creates sorted max file.
 * @author mullermw
 *
 */
public class MergeTask implements Task {
	
	private static Logger logger = Logger.getLogger(MergeTask.class.getName());

	private File outputDir;
	private int[] lengthsOfSplitReads;
	private ReadMask[] readMasks;
	private int[] countsMaxErrorsAllowedWhenMatching;
	private int extensionMatchScore;
	private int extensionMismatchPenalty;
	private File[] filesReferenceFastaPartitions;
	private File[] filesSplitReadsCSFasta;
	private JobSubmissionParameters parameterTemplate;
	private Long maxMemorytPerJobInBytes;
	private Double memoryRequirementAdjustmentFactor;
	private Integer maxMappingLocationsBeforeNotReportingRead;
	private Integer numberOfReadsPerMergeJobOverride;
	private int numberOfReads;
	private Integer minAlignmentScoreForReportingAlignment;
	private Integer uniquenessGap;
	private Integer minMappingLocationsRequiredBeforeReportingRead;
	private Integer lengthOfReads;
	private Integer knownJunctionPenalty;
	private Integer putativeJunctionPenalty;
	private FilteringMode filteringMode;
	private List<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
	private File scratchDir;
	private Boolean deleteIntermediateFiles;
	private Boolean compressIntermediateFiles;
	private File fullLengthReadsCsfasta;
	private File junctionMaxFile;
	
	
	@Override
	/**
	 * Note, the *.max file and *.max.idx file are renamed to *.merged.max
	 * and *.max.idx to indicate that they have been merged. 
	 */
	public void doTask() throws TaskException {
        logger.info("Started merge of MAX files");
        try {
        	
        	FileUtils.deleteDir(getMergeDir());
        	getFileMappingReportTxt().delete();
        	getFileFullyMergedAndSortedMAXcsfasta().delete();
        	getFileFullyMergedAndSortedMAXgff().delete();
        	getFileFullyMergedMAXcsfasta().delete();
        	getFileMappingsHistogramCumulativePDF().delete();
        	getFileMappingsHistogramPDF().delete();
        	getFileReadsWithoutMinScoringAlignmentCSfasta().delete();
        	getMasterJobRemovalScript().delete();
        	getMasterJobSubmissionScript().delete();
        	getMasterListOfJobOutputFiles().delete();
        	
	        ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();
	
	        int maxMinPossibleAlignmentScore = Integer.MIN_VALUE;
	        for (int indexReadSplit = 0; indexReadSplit < lengthsOfSplitReads.length; indexReadSplit++) {
	            int numberOfUnmaskedPositions = readMasks[indexReadSplit].numberOfUnmaskedPositions(lengthsOfSplitReads[indexReadSplit]);
	            maxMinPossibleAlignmentScore = Math.max(maxMinPossibleAlignmentScore,
	                                            (numberOfUnmaskedPositions - countsMaxErrorsAllowedWhenMatching[indexReadSplit]) * extensionMatchScore + countsMaxErrorsAllowedWhenMatching[indexReadSplit] * extensionMismatchPenalty );
	        }
	        logger.info("maxMinPossibleAlignmentScore:\t" + maxMinPossibleAlignmentScore);
	
	        // locate all the MAX files and merge them at the reads file partition level, leaving reference file partitions and read sequence splits intact
	        File[][] filesMAXByRefPartitionByReadSequenceSplit = new File[filesReferenceFastaPartitions.length][filesSplitReadsCSFasta.length];
	        File[][] filesMAXIndexByRefPartitionByReadSequenceSplit = new File[filesReferenceFastaPartitions.length][filesSplitReadsCSFasta.length];
	        for (int indexReferenceFilePartition = 0; indexReferenceFilePartition < filesMAXByRefPartitionByReadSequenceSplit.length; indexReferenceFilePartition++)
	            for (int indexReadsFilePartition = 0; indexReadsFilePartition < filesMAXByRefPartitionByReadSequenceSplit[0].length; indexReadsFilePartition++) {
	
	                File folderCurrentPartition = new File(getMappingAndExtensionDir(), "reads_" + indexReadsFilePartition + "/ref_" + indexReferenceFilePartition + "/");
	                
	                // does the MAX file merged across partitions already exist?
	                File filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[] = folderCurrentPartition.listFiles(WholeTranscriptomeAnalyzer.FILENAME_FILTER_MERGED_MAX);
	                if (filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged.length == 1)
	                    filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[0];
	                else {
	                    // don't know why there would be more than one file, but if so cleanup
	                    for (int i = 0; i < filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged.length; i++)
	                        filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[i].delete();
	
	                    // if not, has it been compressed?
	                    File filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[] = folderCurrentPartition.listFiles(WholeTranscriptomeAnalyzer.FILENAME_FILTER_MERGED_MAX_COMPRESSED);
	                    if (filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped.length == 1) {
	                        logger.info("Decompressing intermediate files...");
	                        CompressionUtilities.decompressFiles(filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[0], folderCurrentPartition);
	                        filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = folderCurrentPartition.listFiles(WholeTranscriptomeAnalyzer.FILENAME_FILTER_MERGED_MAX)[0];
	                        logger.info("OK");
	                    } else {
	
	                        // don't know why there would be more than one file, but if so cleanup
	                        for (int i = 0; i < filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped.length; i++)
	                            filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[i].delete();
	
	                        // if not compressed, are the unmerged MAX files available to make the merge anew?
	                        File filesMAXForThisRefPartitionAndThisReadSequenceSplit[] = folderCurrentPartition.listFiles(WholeTranscriptomeAnalyzer.FILENAME_FILTER_MAX);
	                        if (filesMAXForThisRefPartitionAndThisReadSequenceSplit.length == 1) {
	                            int indexOfSecondToLastPeriod = filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().lastIndexOf(".", filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().length() - WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE.length() -1);
	                            filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().substring(0, indexOfSecondToLastPeriod) + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
	                            filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition].getPath() + WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE);
	                            filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].renameTo(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);
	                            new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath() + WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE).renameTo(filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);
	
	                        } else if (filesMAXForThisRefPartitionAndThisReadSequenceSplit.length > 1) {
	
	                            throw new Exception("MERGE of multiple reads file partitioned MAX files is not currently supported.");
	
	                            /*System.out.println(" -->Merging MAX files, by reads file partition:\t" + folderCurrentPartition);
	                            int indexOfSecondToLastPeriod = filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().lastIndexOf(".", filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().length() - EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE.length() -1);
	                            filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().substring(0, indexOfSecondToLastPeriod) + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
	                            filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition].getPath() + EXTENSION_FOR_INDEX_FILE);
	                            Arrays.sort(filesMAXForThisRefPartitionAndThisReadSequenceSplit, new com.lifetechnologies.solid.wt.MAFileNameComparator());
	                            com.lifetechnologies.solid.wt.Utilities.concatenateFiles(filesMAXForThisRefPartitionAndThisReadSequenceSplit, filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition], true);
	
	                            File filesMAXIndexForThisRefPartitionAndThisReadSequenceSplit[] = new File[filesMAXForThisRefPartitionAndThisReadSequenceSplit.length];
	                            for (int i = 0; i < filesMAXIndexForThisRefPartitionAndThisReadSequenceSplit.length; i++)
	                                filesMAXIndexForThisRefPartitionAndThisReadSequenceSplit[i] = new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[i].getPath() + EXTENSION_FOR_INDEX_FILE);
	                            // FIXME: This is not the proper way to merge index files and is a temporary hack that works only when there is a single readsFilePartition 
	                            com.lifetechnologies.solid.wt.Utilities.concatenateFiles(filesMAXIndexForThisRefPartitionAndThisReadSequenceSplit, filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition], true);
	
	                            for (int i = 0; i < filesMAXForThisRefPartitionAndThisReadSequenceSplit.length; i++)
	                                filesMAXForThisRefPartitionAndThisReadSequenceSplit[i].delete();
	                              */
	                        } else
	                            throw new Exception("Could not find MAX files in:\t" + folderCurrentPartition.getPath() + "\nPlease re-run with RUN_EXTENSION turned on.");
	
	                    }
	                }
	                listOfFilesThatAreNoLongerNeeded.add(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);
	            }
	
	        JobSubmissionParameters params = parameterTemplate.clone();
	        params.setMemoryRequirement(maxMemorytPerJobInBytes.longValue());
	        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
	        // Need to break up the reads into a new set of partitions that can fit into memory on each node (for sorting).  This is a worst case guess for now.
	        // FIXME: IS the calculation robust?
	        int numberOfReadsPerPartition = (int) (maxMemorytPerJobInBytes/memoryRequirementAdjustmentFactor / 4000.0 *  10.0 / maxMappingLocationsBeforeNotReportingRead);   //1500.0);
	        if (numberOfReadsPerMergeJobOverride != null)
	            numberOfReadsPerPartition = numberOfReadsPerMergeJobOverride;
	        File filesOutputForEachReadsFilePartition[] = new File[(int)Math.ceil((double)numberOfReads / numberOfReadsPerPartition)];
	        File filesSortedOutputForEachReadsFilePartition[] = new File[filesOutputForEachReadsFilePartition.length];
	        File filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[] = new File[filesOutputForEachReadsFilePartition.length];
	        File filesSerializedReportOfExtendedReadMappingsOut[] = new File[filesOutputForEachReadsFilePartition.length];
	        for (int indexOfPartition = 0; indexOfPartition < filesOutputForEachReadsFilePartition.length; indexOfPartition++) {
	
	            int indexOfFirstReadToMerge = indexOfPartition * numberOfReadsPerPartition;
	            int indexOfLastReadToMerge = Math.min(indexOfFirstReadToMerge + numberOfReadsPerPartition -1, numberOfReads -1);
	
	            logger.info("-->Merging MAX files by reference partition and read sequence split for reads " + indexOfFirstReadToMerge + " to " + indexOfLastReadToMerge);
	
	            File folderForMergingThisReadPartition = new File(getMergeDir(), "reads_" + indexOfFirstReadToMerge + "_to_" + indexOfLastReadToMerge);
	            if (!folderForMergingThisReadPartition.exists())
	                folderForMergingThisReadPartition.mkdirs();
	
	            filesOutputForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE.substring(1) + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	            listOfFilesThatAreNoLongerNeeded.add(filesOutputForEachReadsFilePartition[indexOfPartition]);
	
	            filesSortedOutputForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + WholeTranscriptomeAnalyzer.EXTENSION_FOR_SORTED_OUTPUT_FILE.substring(1) + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	            listOfFilesThatAreNoLongerNeeded.add(filesSortedOutputForEachReadsFilePartition[indexOfPartition]);
	
	            filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_WITHOUT_MIN_SCORING_ALIGNMENT_FILE.substring(1) + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	            listOfFilesThatAreNoLongerNeeded.add(filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition]);
	
	            filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/reportOfExtendedReadMappingsOut." + System.currentTimeMillis() + ".obj");
	            listOfFilesThatAreNoLongerNeeded.add(filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition]);
	
	            ExtendedReadMappingsMerger.startMergeExtendedReadMappingsJob(minAlignmentScoreForReportingAlignment,
	                                                                        uniquenessGap,
	                                                                        minMappingLocationsRequiredBeforeReportingRead,
	                                                                        maxMappingLocationsBeforeNotReportingRead,
	                                                                        maxMinPossibleAlignmentScore,
	                                                                        lengthOfReads * extensionMatchScore,
	                                                                        indexOfFirstReadToMerge,
	                                                                        indexOfLastReadToMerge,
	                                                                        knownJunctionPenalty,
	                                                                        putativeJunctionPenalty,
	                                                                        this.filteringMode, 
	                                                                        WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE,
	                                                                        listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
	                                                                        filesMAXByRefPartitionByReadSequenceSplit,
	                                                                        this.junctionMaxFile,
	                                                                        //reportOfExtendedReadMappings,
	                                                                        clusterInterface,
	                                                                        "wt_merge.",
	                                                                        WholeTranscriptomeAnalyzer.getLibDir().getPath() + "/'*':"+WholeTranscriptomeAnalyzer.getWholeTranscriptomeJar().getPath(),
	                                                                        scratchDir,
	                                                                        folderForMergingThisReadPartition,
	                                                                        folderForMergingThisReadPartition,
	                                                                        filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition],
	                                                                        filesOutputForEachReadsFilePartition[indexOfPartition],
	                                                                        filesSortedOutputForEachReadsFilePartition[indexOfPartition],
	                                                                        filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition],
	                                                                        true);
	
	
	        }
	
	        File fileMasterMergeJobSubmissionScript = getMasterJobSubmissionScript();
	        clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterMergeJobSubmissionScript);
	        File fileMasterMergeJobRemovalScript = getMasterJobRemovalScript();
	        clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterMergeJobRemovalScript);
	        File fileMasterMergeListOfScriptOutputFiles = getMasterListOfJobOutputFiles();
	        clusterInterface.writeMasterListOfJobOutputFilesFromLog(fileMasterMergeListOfScriptOutputFiles);
	
	        fileMasterMergeJobSubmissionScript.setExecutable(true, true);
	        fileMasterMergeJobRemovalScript.setExecutable(true, true);
	
	        // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
	        logger.info("Waiting for merge jobs to finish...");
	        while (!clusterInterface.checkIfLoggedJobsComplete())
	            Thread.sleep(5000);
	
	        logger.info("OK");
	
	        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
	        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
	            StringBuffer msg = new StringBuffer("The following jobs failed:\n");
	            Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
	            while (iteratorOverFailFiles.hasNext()) {
	                msg.append(( iteratorOverFailFiles.next()).getPath());
	                msg.append("\n");
	            }
	            logger.info(msg.toString());
	            throw new Exception("One or more merge jobs failed.");
	        }
	        
	
	        logger.info(" -->Merging unsorted MAX files, by reads file partition.");
	        Utilities.concatenateFiles(filesOutputForEachReadsFilePartition, getFileFullyMergedMAXcsfasta(), false);
	        ExtendedReadMappingsReport reportOfExtendedReadMappings = ExtendedReadMappingsReport.combineSerializedMappingReports(filesSerializedReportOfExtendedReadMappingsOut);
	        reportOfExtendedReadMappings.setFilteringMode(filteringMode);
	
	        logger.info(" -->Merging sorted MAX files, by reads file partition.");
	        BufferedWriter writer = new BufferedWriter(new FileWriter(getFileFullyMergedAndSortedMAXcsfasta()));
	        BufferedReader readersSortedOutputFilesByReadsFilePartition[] = new BufferedReader[filesSortedOutputForEachReadsFilePartition.length];
	        for (int i = 0; i < filesSortedOutputForEachReadsFilePartition.length; i++)
	            readersSortedOutputFilesByReadsFilePartition[i] = new BufferedReader(new FileReader(filesSortedOutputForEachReadsFilePartition[i]));
	
	        TreeMap<ExtendedReadMapping, Integer> mapExtendedReadMappingToPartitionIndex = new TreeMap<ExtendedReadMapping, Integer>(new ReadMappingComparator());
	        String lines[] = new String[filesOutputForEachReadsFilePartition.length];
	        for (int indexReadsFilePartition = 0; indexReadsFilePartition < filesSortedOutputForEachReadsFilePartition.length; indexReadsFilePartition++) {
	            lines[indexReadsFilePartition] = readersSortedOutputFilesByReadsFilePartition[indexReadsFilePartition].readLine();
	            if (lines[indexReadsFilePartition] != null) {
	                ExtendedReadMapping extendedReadMappings[] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(lines[indexReadsFilePartition].substring(1));
	                mapExtendedReadMappingToPartitionIndex.put(extendedReadMappings[0], indexReadsFilePartition);
	            }
	        }
	
	        while (mapExtendedReadMappingToPartitionIndex.size() > 0) {
	            ExtendedReadMapping extendedReadMappingFirst = mapExtendedReadMappingToPartitionIndex.firstKey();
	            Integer indexOfCorrespondingPartition = mapExtendedReadMappingToPartitionIndex.get(extendedReadMappingFirst);
	            mapExtendedReadMappingToPartitionIndex.remove(extendedReadMappingFirst);
	            writer.write(">" + extendedReadMappingFirst.getIdOfRead() + extendedReadMappingFirst.toString());
	            writer.newLine();
	            writer.write(readersSortedOutputFilesByReadsFilePartition[indexOfCorrespondingPartition].readLine());
	            writer.newLine();
	            lines[indexOfCorrespondingPartition] = readersSortedOutputFilesByReadsFilePartition[indexOfCorrespondingPartition].readLine();
	            if (lines[indexOfCorrespondingPartition] != null) {
	                ExtendedReadMapping extendedReadMappings[] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(lines[indexOfCorrespondingPartition].substring(1));
	                mapExtendedReadMappingToPartitionIndex.put(extendedReadMappings[0], indexOfCorrespondingPartition);
	            }
	        }
	
	        for (int i = 0; i < filesSortedOutputForEachReadsFilePartition.length; i++)
	            readersSortedOutputFilesByReadsFilePartition[i].close();
	        writer.close();
	        
	        logger.info("Merging sorted GFF files, by reads file partition.");
	        Collection<Reader> gffSrc = new ArrayList<Reader>();
	        PrintWriter gffWriter = null;
	        try {
	        	for (File file : filesSortedOutputForEachReadsFilePartition ) 
	        		gffSrc.add(new BufferedReader(new FileReader(new File(file.getParent(), file.getName().replace(".csfasta", ".gff")))));
	        	Date now = new Date();
	        	gffWriter = new PrintWriter(new BufferedWriter(new FileWriter(getFileFullyMergedAndSortedMAXgff())));
	        	gffWriter.println("##gff-version 2");
	        	gffWriter.printf ("##source-version WTP v%s\n", WTConfigUtils.getApplicationProperties().getString("version"));
	        	gffWriter.printf ("##date %tY-%tb-%td\n", now, now, now);
	        	gffWriter.printf ("##time %tH:%tM:%tS\n", now, now, now);
	        	gffWriter.println("##Type solid_read");
	        	gffWriter.println("##color-code AA=0,AC=1,AG=2,AT=3,CA=1,CC=0,CG=3,CT=2,GA=2,GC=3,GG=0,GT=1,TA=3,TC=2,TG=1,TT=0");
	        	gffWriter.println("##primer-base T");
	        	gffWriter.println("##attribute definition");
	        	gffWriter.println("##bd bead id");
	        	gffWriter.println("##rs position of alignment start in read, starting with 1.");
	        	gffWriter.println("##mm number of mismatches in the alignment.");
	        	gffWriter.println("##g color sequence of read.");
	        	gffWriter.println("##i index of reference sequence for this alignment, starting with 1.");
	        	gffWriter.println("##jp presence indicates this read maps to a splice junction.  Value is the start");
	        	gffWriter.println("##   position of the other feature completing this alignment.");
	        	gffWriter.println("##jt type of splice junction aligned to, k=known, p=putative");
		        LineUtils.mergeLines(gffSrc, gffWriter, GffFeature.FORMAT);
	        } finally {
	        	for (Reader reader : gffSrc) reader.close();
	        	if (gffWriter != null) gffWriter.close();
	        }
	
	        logger.info(" -->Merging reads w/o min scoring alignments files, by reads file partition.");
	        Utilities.concatenateFiles(filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition, getFileReadsWithoutMinScoringAlignmentCSfasta(), false);
	
	
	        if (deleteIntermediateFiles || compressIntermediateFiles) {
	            if (compressIntermediateFiles && deleteIntermediateFiles)  logger.info("Compressing and deleting intermediate files...");
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
	            logger.info("OK");
	        }
	
	
	       logger.info("\nFinished merge of MAX files");
	
	        try {
	            reportOfExtendedReadMappings.getHistogramOfMappingsByScore().setReportRelativeFrequencies(true);
	            reportOfExtendedReadMappings.getHistogramOfMappingsByScore().toPDF(getFileMappingsHistogramPDF(), false, false, false, null);
	            reportOfExtendedReadMappings.getHistogramOfMappingsByScore().toPDF(getFileMappingsHistogramCumulativePDF(), false, true, false, null);
	
	        } catch (Exception e) {
	            logger.log(Level.SEVERE, "Caught exception while trying to generate PDF.", e);
	            e.printStackTrace();
	        }
	
	        reportOfExtendedReadMappings.writeReportToTextFile(getFileMappingReportTxt());
	        logger.info("\n"+reportOfExtendedReadMappings);
	        
        } catch (Exception e) {
        	logger.log(Level.SEVERE, e.getMessage(), e);
        	throw new TaskException(e);
        }
        
    }


	@Override
	public boolean isDone() throws TaskException {
		// TODO Auto-generated method stub
		return false;
	}
	
	@Override
	public boolean isRemote() {
		return true;
	}
	
	@Override
	public void setRemote(boolean remote) {
		if (remote != true) throw new IllegalArgumentException("MergeTask only supported as remote task.");
	}

	public File getOutputDir() {
		return outputDir;
	}

	public void setOutputDir(File outputDir) {
		this.outputDir = outputDir;
	}

	public int[] getLengthsOfSplitReads() {
		return lengthsOfSplitReads;
	}


	public void setLengthsOfSplitReads(int[] lengthsOfSplitReads) {
		this.lengthsOfSplitReads = lengthsOfSplitReads;
	}


	public ReadMask[] getReadMasks() {
		return readMasks;
	}


	public void setReadMasks(ReadMask[] readMasks) {
		this.readMasks = readMasks;
	}


	public int[] getCountsMaxErrorsAllowedWhenMatching() {
		return countsMaxErrorsAllowedWhenMatching;
	}


	public void setCountsMaxErrorsAllowedWhenMatching(
			int[] countsMaxErrorsAllowedWhenMatching) {
		this.countsMaxErrorsAllowedWhenMatching = countsMaxErrorsAllowedWhenMatching;
	}


	public int getExtensionMatchScore() {
		return extensionMatchScore;
	}


	public void setExtensionMatchScore(int extensionMatchScore) {
		this.extensionMatchScore = extensionMatchScore;
	}


	public int getExtensionMismatchPenalty() {
		return extensionMismatchPenalty;
	}


	public void setExtensionMismatchPenalty(int extensionMismatchPenalty) {
		this.extensionMismatchPenalty = extensionMismatchPenalty;
	}


	public File[] getFilesReferenceFastaPartitions() {
		return filesReferenceFastaPartitions;
	}


	public void setFilesReferenceFastaPartitions(
			File[] filesReferenceFastaPartitions) {
		this.filesReferenceFastaPartitions = filesReferenceFastaPartitions;
	}


	public File[] getFilesSplitReadsCSFasta() {
		return filesSplitReadsCSFasta;
	}


	public void setFilesSplitReadsCSFasta(File[] filesSplitReadsCSFasta) {
		this.filesSplitReadsCSFasta = filesSplitReadsCSFasta;
	}


	public JobSubmissionParameters getParameterTemplate() {
		return parameterTemplate;
	}


	public void setParameterTemplate(JobSubmissionParameters parameterTemplate) {
		this.parameterTemplate = parameterTemplate;
	}


	public Long getMaxMemorytPerJobInBytes() {
		return maxMemorytPerJobInBytes;
	}


	public void setMaxMemorytPerJobInBytes(Long maxMemorytPerJobInBytes) {
		this.maxMemorytPerJobInBytes = maxMemorytPerJobInBytes;
	}


	public Double getMemoryRequirementAdjustmentFactor() {
		return memoryRequirementAdjustmentFactor;
	}


	public void setMemoryRequirementAdjustmentFactor(
			Double memoryRequirementAdjustmentFactor) {
		this.memoryRequirementAdjustmentFactor = memoryRequirementAdjustmentFactor;
	}


	public Integer getMaxMappingLocationsBeforeNotReportingRead() {
		return maxMappingLocationsBeforeNotReportingRead;
	}


	public void setMaxMappingLocationsBeforeNotReportingRead(
			Integer maxMappingLocationsBeforeNotReportingRead) {
		this.maxMappingLocationsBeforeNotReportingRead = maxMappingLocationsBeforeNotReportingRead;
	}


	public Integer getNumberOfReadsPerMergeJobOverride() {
		return numberOfReadsPerMergeJobOverride;
	}


	public void setNumberOfReadsPerMergeJobOverride(
			Integer numberOfReadsPerMergeJobOverride) {
		this.numberOfReadsPerMergeJobOverride = numberOfReadsPerMergeJobOverride;
	}


	public int getNumberOfReads() {
		return numberOfReads;
	}


	public void setNumberOfReads(int numberOfReads) {
		this.numberOfReads = numberOfReads;
	}


	public Integer getMinAlignmentScoreForReportingAlignment() {
		return minAlignmentScoreForReportingAlignment;
	}


	public void setMinAlignmentScoreForReportingAlignment(
			Integer minAlignmentScoreForReportingAlignment) {
		this.minAlignmentScoreForReportingAlignment = minAlignmentScoreForReportingAlignment;
	}


	public Integer getUniquenessGap() {
		return uniquenessGap;
	}


	public void setUniquenessGap(Integer uniquenessGap) {
		this.uniquenessGap = uniquenessGap;
	}


	public Integer getMinMappingLocationsRequiredBeforeReportingRead() {
		return minMappingLocationsRequiredBeforeReportingRead;
	}


	public void setMinMappingLocationsRequiredBeforeReportingRead(
			Integer minMappingLocationsRequiredBeforeReportingRead) {
		this.minMappingLocationsRequiredBeforeReportingRead = minMappingLocationsRequiredBeforeReportingRead;
	}


	public Integer getLengthOfReads() {
		return lengthOfReads;
	}


	public void setLengthOfReads(Integer lengthOfReads) {
		this.lengthOfReads = lengthOfReads;
	}


	public Integer getKnownJunctionPenalty() {
		return knownJunctionPenalty;
	}


	public void setKnownJunctionPenalty(Integer knownJunctionPenalty) {
		this.knownJunctionPenalty = knownJunctionPenalty;
	}


	public Integer getPutativeJunctionPenalty() {
		return putativeJunctionPenalty;
	}


	public void setPutativeJunctionPenalty(Integer putativeJunctionPenalty) {
		this.putativeJunctionPenalty = putativeJunctionPenalty;
	}


	public FilteringMode getFilteringMode() {
		return filteringMode;
	}


	public void setFilteringMode(FilteringMode filteringMode) {
		this.filteringMode = filteringMode;
	}


	public List<Integer>[] getListsOfFullRefSequenceNumbersOrderedForEachReferencePartition() {
		return listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
	}


	public void setListsOfFullRefSequenceNumbersOrderedForEachReferencePartition(
			List<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition) {
		this.listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
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


	public File getFullLengthReadsCsfasta() {
		return fullLengthReadsCsfasta;
	}


	public void setFullLengthReadsCsfasta(File fullLengthReadsCsfasta) {
		this.fullLengthReadsCsfasta = fullLengthReadsCsfasta;
	}


	public File getJunctionMaxFile() {
		return junctionMaxFile;
	}


	public void setJunctionMaxFile(File junctionMaxFile) {
		this.junctionMaxFile = junctionMaxFile;
	}


	public File getTmpDir() {
		return new File(this.getOutputDir(), "tmp");
	}
	
	public File getScriptDir() {
		return new File(this.getOutputDir(), "scripts");
	}
	
	public File getResultDir() {
		return new File(this.getOutputDir(), "output");
	}
	
	public File getMappingAndExtensionDir() {
		return new File(this.getTmpDir(), "mappingAndExtension");
	}
	
	public File getMergeDir() {
		return new File(this.getTmpDir(), "merge");
	}
	
	public File getMasterJobSubmissionScript() {
		return new File(getScriptDir(), "/merge.master.script_submission.sh");
	}
	
	public File getMasterJobRemovalScript() {
		return new File(getScriptDir(), "merge.master.script_removal.sh");
	}
	
	public File getMasterListOfJobOutputFiles() {
		return new File(getScriptDir(), "merge.list_of_script_output_files.out");
	}
	
	public File getFileFullyMergedMAXcsfasta() {
		return new File(getResultDir(), fullLengthReadsCsfasta.getName().substring(0, fullLengthReadsCsfasta.getName().lastIndexOf('.'))
                + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	}
	
	public File getFileFullyMergedAndSortedMAXcsfasta() {
		return new File(getResultDir(), fullLengthReadsCsfasta.getName().substring(0, fullLengthReadsCsfasta.getName().lastIndexOf('.'))
                + WholeTranscriptomeAnalyzer.EXTENSION_FOR_SORTED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	}
	
	public File getFileFullyMergedAndSortedMAXgff() {
		return new File(getResultDir(), fullLengthReadsCsfasta.getName().substring(0, fullLengthReadsCsfasta.getName().lastIndexOf('.'))
                + WholeTranscriptomeAnalyzer.EXTENSION_FOR_SORTED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_MERGED_OUTPUT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_FILTERED_OUTPUT_FILE + ".gff");
	}
	
	public File getFileReadsWithoutMinScoringAlignmentCSfasta() {
		return new File(getResultDir(), fullLengthReadsCsfasta.getName().substring(0, fullLengthReadsCsfasta.getName().lastIndexOf('.'))
                + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_WITHOUT_MIN_SCORING_ALIGNMENT_FILE + WholeTranscriptomeAnalyzer.EXTENSION_FOR_READS_FILE);
	}
	
	public File getFileMappingReportTxt() {
		return new File(getResultDir(), "alignmentReport.txt");
	}
	
	public File getFileMappingsHistogramPDF() {
		return new File(getResultDir(), "alignmentsByScore.histograms.pdf");
	}
	
	public File getFileMappingsHistogramCumulativePDF() {
		return new File(getResultDir(), "alignmentsByScore.histograms.cumulative.pdf");
	}

	
}
