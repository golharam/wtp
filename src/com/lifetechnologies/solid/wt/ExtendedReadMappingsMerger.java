package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.mapper.FilteringMode;
import com.lifetechnologies.solid.wt.mapper.GffFeature;
import com.lifetechnologies.solid.wt.splice.FlankedJunction;
import com.lifetechnologies.solid.wt.splice.JunctionFastaFile;
import com.lifetechnologies.util.LineUtils;


/**
 * User: tuchbb
 * Date: Oct 15, 2008
 * Time: 8:55:19 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 * 
 */
public class ExtendedReadMappingsMerger {

	static final Logger logger = Logger.getLogger(ExtendedReadMappingsMerger.class.toString());
    public static String PATH_TO_JAVA_EXE = "java";

    //public static int BUFFER_SIZE_FOR_READING_AND_WRITING_FILES_IN_BYTES = 1000000;
    
    public static void main(String args[]) throws Exception {
    	//TODO remove test code
    	if (System.getProperty("os.name").toLowerCase().contains("windows")) {
    		JunctionMapping.main(args);
    		System.exit(0);
    	}
        int minAlignmentScoreForReportingAlignment = Integer.parseInt(args[0]);
        int minScoreGapToSecondBestAlignmentForUniqueness = Integer.parseInt(args[1]);
        int minMappingsRequiredBeforeReportingRead = Integer.parseInt(args[2]);
        int maxMappingsAllowedBeforeNotReportingRead = Integer.parseInt(args[3]);
        int minPossibleAlignmentScore = Integer.parseInt(args[4]);
        int maxPossibleAlignmentScore = Integer.parseInt(args[5]);
        int indexOfFirstReadToMerge = Integer.parseInt(args[6]);
        int indexOfLastReadToMerge = Integer.parseInt(args[7]);
        int knownJunctionPenalty = Integer.valueOf(args[8]);
        int putativeJunctionPenalty = Integer.valueOf(args[9]);

        boolean filterReadsWithOneOrMoreSplitsHavingAFilterTag = Boolean.parseBoolean(args[10]);

        String extensionForMAXFileIndexFiles = args[11];

        ObjectInputStream objectInputStream = new ObjectInputStream(new FileInputStream(args[12]));
        @SuppressWarnings("unchecked")
        ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = (ArrayList<Integer>[])objectInputStream.readObject();
        objectInputStream.close();

        objectInputStream = new ObjectInputStream(new FileInputStream(args[13]));
        File[][] filesMAXByRefPartitionByReadSequenceSplit = (File[][])objectInputStream.readObject();
        objectInputStream.close();

        objectInputStream = new ObjectInputStream(new FileInputStream(args[14]));
        ExtendedReadMappingsReport reportOfExtendedReadMappings = (ExtendedReadMappingsReport)objectInputStream.readObject();
        objectInputStream.close();

        String outputReportFile = args[15];
        
        File fileMergedOutput = new File(args[16]);
        File fileMergedSortedOutput = new File(args[17]);
        File fileReadsWithoutMinScoringAlignment = new File(args[18]);
        
        File fileJunctionMax = null;
        if (args.length > 19) fileJunctionMax = new File(args[19]);
        
        mergeMAXFilePartitionsAndReadSequenceSplits(minAlignmentScoreForReportingAlignment, minScoreGapToSecondBestAlignmentForUniqueness,
                                                    minMappingsRequiredBeforeReportingRead, maxMappingsAllowedBeforeNotReportingRead,
                                                    minPossibleAlignmentScore, maxPossibleAlignmentScore,
                                                    indexOfFirstReadToMerge, indexOfLastReadToMerge,
                                                    knownJunctionPenalty, putativeJunctionPenalty,
                                                    filterReadsWithOneOrMoreSplitsHavingAFilterTag, extensionForMAXFileIndexFiles,
                                                    listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
                                                    filesMAXByRefPartitionByReadSequenceSplit, fileJunctionMax,
                                                    reportOfExtendedReadMappings, fileMergedOutput, fileMergedSortedOutput, fileReadsWithoutMinScoringAlignment);

        ObjectOutputStream objectOutputStream = new ObjectOutputStream(new FileOutputStream(outputReportFile));
        objectOutputStream.writeObject(reportOfExtendedReadMappings);
        objectOutputStream.flush();
        objectOutputStream.close();
    }


    public static File startMergeExtendedReadMappingsJob(int minAlignmentScoreForReportingAlignment,
                                                         int minScoreGapToSecondBestAlignmentForUniqueness,
                                                         int minMappingsRequiredBeforeReportingRead,
                                                         int maxMappingsAllowedBeforeNotReportingRead,
                                                         int minPossibleAlignmentScore,
                                                         int maxPossibleAlignmentScore,
                                                         int indexOfFirstReadToMerge,
                                                         int indexOfLastReadToMerge,
                                                         int knownJunctionPenalty,
                                                         int putativeJunctionPenalty,
                                                         FilteringMode filteringMode,
                                                         String extensionForMAXFileIndexFiles,
                                                         List<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
                                                         File[][] filesMAXByRefPartitionByReadSequenceSplit,
                                                         File fileJunctionMAX,
                                                         ClusterInterface clusterInterface,
                                                         String prefixForJobNames,
                                                         String classPathForJavaCode,
                                                         File folderForTempFilesOnComputeNodes,
                                                         File folderForScriptIO,
                                                         File folderForSerializedObjects,
                                                         File fileSerializedReportOfExtendedReadMappingsOut,
                                                         File fileMergedOutput,
                                                         File fileMergedSortedOutput,
                                                         File fileReadsWithoutMinScoringAlignment,
                                                         boolean verbose) throws Exception {

    	StringBuffer buffer = new StringBuffer();
    	buffer.append("\n");
        buffer.append("---------------------------------------------------------------------------------------------\n");
        buffer.append("merged output file:\t" + fileMergedOutput.getPath()+"\n");
        buffer.append("min alignment score for reporting alignment:\t" + minAlignmentScoreForReportingAlignment+"\n");
        buffer.append("min mappings required before reporting read:\t" + minMappingsRequiredBeforeReportingRead+"\n");
        buffer.append("max mappings allowed before not reporting read:\t" + maxMappingsAllowedBeforeNotReportingRead+"\n");
        buffer.append("index of first read to merge:\t" + indexOfFirstReadToMerge+"\n");
        buffer.append("index of last read to merge:\t" + indexOfLastReadToMerge+"\n");
        buffer.append("---------------------------------------------------------------------------------------------\n");
        logger.info(buffer.toString());
        
        ArrayList<String> commandStrings = new ArrayList<String>();

        File fileSerializedListsOfFullRefSequenceNumbersOrderedForEachReferencePartition = new File(folderForSerializedObjects.getAbsolutePath() + "/mapOfSequenceNumbers." + System.currentTimeMillis() + ".obj");
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(new FileOutputStream(fileSerializedListsOfFullRefSequenceNumbersOrderedForEachReferencePartition));
        objectOutputStream.writeObject(listsOfFullRefSequenceNumbersOrderedForEachReferencePartition);
        objectOutputStream.flush();
        objectOutputStream.close();

        //Convert File paths to absolute paths.
        File[][] arrCopy = new File[filesMAXByRefPartitionByReadSequenceSplit.length][];
        for (int i = 0; i<arrCopy.length; i++ ) {
        	arrCopy[i] = new File[filesMAXByRefPartitionByReadSequenceSplit[i].length];
        	for (int j=0; j<arrCopy[i].length; j++) {
        		arrCopy[i][j] = filesMAXByRefPartitionByReadSequenceSplit[i][j].getAbsoluteFile();
        	}
        }
        filesMAXByRefPartitionByReadSequenceSplit = arrCopy;
        File fileSerializedFilesMAXByRefPartitionByReadSequenceSplit = new File(folderForSerializedObjects.getAbsolutePath() + "/arrayOfMaxFiles." + System.currentTimeMillis() + ".obj");
        objectOutputStream = new ObjectOutputStream(new FileOutputStream(fileSerializedFilesMAXByRefPartitionByReadSequenceSplit));
        objectOutputStream.writeObject(filesMAXByRefPartitionByReadSequenceSplit);
        objectOutputStream.flush();
        objectOutputStream.close();

        File fileSerializedReportOfExtendedReadMappingsIn = new File(folderForSerializedObjects.getAbsolutePath() + "/reportOfExtendedReadMappingsIn." + System.currentTimeMillis() + ".obj");
        objectOutputStream = new ObjectOutputStream(new FileOutputStream(fileSerializedReportOfExtendedReadMappingsIn));
        objectOutputStream.writeObject(new ExtendedReadMappingsReport());
        objectOutputStream.flush();
        objectOutputStream.close();

        String command = PATH_TO_JAVA_EXE + 
        					" -classpath " +  classPathForJavaCode + 
        					" -Xmx" + (int)Math.ceil(clusterInterface.getJobSubmissionParameters().getMemoryRequirement() / Constants.BYTES_PER_MEGABYTE) + "m" +
        					(Utilities.isTrue(System.getProperty(WholeTranscriptomeAnalyzer.EXPERIMENTAL_MODE_SYSTEM_PROPERTY)) ? " -D"+WholeTranscriptomeAnalyzer.EXPERIMENTAL_MODE_SYSTEM_PROPERTY+"=true" : "") +
        					" " + ExtendedReadMappingsMerger.class.getName() +
                            " " + minAlignmentScoreForReportingAlignment +
                            " " + minScoreGapToSecondBestAlignmentForUniqueness +
                            " " + minMappingsRequiredBeforeReportingRead +
                            " " + maxMappingsAllowedBeforeNotReportingRead +
                            " " + minPossibleAlignmentScore +
                            " " + maxPossibleAlignmentScore +
                            " " + indexOfFirstReadToMerge +
                            " " + indexOfLastReadToMerge +
                            " " + knownJunctionPenalty +
                            " " + putativeJunctionPenalty +
                            " " + (filteringMode == FilteringMode.ONE_OR_MORE) +
                            " " + extensionForMAXFileIndexFiles + 
                            " " + fileSerializedListsOfFullRefSequenceNumbersOrderedForEachReferencePartition.getAbsolutePath() +
                            " " + fileSerializedFilesMAXByRefPartitionByReadSequenceSplit.getAbsolutePath() +
                            " " + fileSerializedReportOfExtendedReadMappingsIn.getAbsolutePath() +
                            " " + fileSerializedReportOfExtendedReadMappingsOut.getAbsolutePath() +
                            " " + fileMergedOutput.getAbsolutePath() +
                            " " + fileMergedSortedOutput.getAbsolutePath() +
                            " " + fileReadsWithoutMinScoringAlignment.getAbsolutePath();

        if (fileJunctionMAX != null) 
        	command += " " + fileJunctionMAX.getAbsolutePath();
        
        

        commandStrings.add(command);

        File fileScript = new File(folderForScriptIO.getPath() + "/" + prefixForJobNames + fileMergedOutput.getName() + ".sh");

        //String workingFolderName = prefixForJobNames.concat(ProcessId.getProcessId())+"."+System.currentTimeMillis();
        //File workingFolder = new File(folderForTempFilesOnComputeNodes, workingFolderName);
        File workingFolder = Utilities.getUniqueFileName(folderForTempFilesOnComputeNodes, prefixForJobNames);
        
        clusterInterface.setLocalTemporaryWorkingPath(workingFolder.getAbsolutePath());

        return clusterInterface.executeJob(fileScript, commandStrings);
    }


    /**
     * Will merge MAX files that have been paritioned by reference and/or by reads, and will handle the merging of individually split read sequences.
     *
     * @param minAlignmentScoreForReportingAlignment
     * @param minMappingsRequiredBeforeReportingRead
     * @param maxMappingsAllowedBeforeNotReportingRead
     * @param maxPossibleAlignmentScore for sizing the histogram
     * @param filterReadsWithOneOrMoreSplitsHavingAFilterTag
     *@param listsOfFullRefSequenceNumbersOrderedForEachReferencePartition  this data structure specifies the following mapping:  reference_partition_index -> sequence_number_in_ref_partition -> sequence_number_in_full_reference
     *                                                                       if not specified, then sequence numbers will be the same in the merged file as in the partitioned MA files (not recommended unless # of ref file paritions == 1).
     * @param filesMAXByRefPartitionByReadSequenceSplit
     * @param fileMergedOutput @throws Exception
     */
    public static void mergeMAXFilePartitionsAndReadSequenceSplits(int minAlignmentScoreForReportingAlignment,
                                                                   int minScoreGapToSecondBestAlignmentForUniqueness,
                                                                    int minMappingsRequiredBeforeReportingRead,
                                                                    int maxMappingsAllowedBeforeNotReportingRead,
                                                                    int minPossibleAlignmentScore,
                                                                    int maxPossibleAlignmentScore,
                                                                    int indexOfFirstReadToMerge,
                                                                    int indexOfLastReadToMerge,
                                                                    int penaltyForKnownJunctions,
                                                                    int penaltyForPutativeJunctions,
                                                                    boolean filterReadsWithOneOrMoreSplitsHavingAFilterTag,
                                                                    String extensionForMAXFileIndexFiles,
                                                                    ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
                                                                    File[][] filesMAXByRefPartitionByReadSequenceSplit,
                                                                    File fileJunctionMAX,
                                                                    ExtendedReadMappingsReport reportOfExtendedReadMappings,
                                                                    File fileMergedOutput,
                                                                    File fileMergedSortedOutput,
                                                                    File fileReadsWithoutMinScoringAlignment) throws Exception {

    	
    	//Capture distributions of scores...will end up in PDFs.
        Histogram histogramOfMappingsByScore = new Histogram(maxPossibleAlignmentScore +1, 1, 0, "Alignment Score");
        String nameOfBestMappingsHistogramSeries = "Best Read Alignments";
        String nameOfUniqueMappingsHistogramSeries = "Unique Read Alignments";
        String nameOfAllMappingsHistogramSeries = "All Read Alignments";
        String nameOfUniqueInclJunctionsHistogramSeries = "Unique Read Alignments including junction alignments";
        histogramOfMappingsByScore.createNewHistogramSeries(nameOfBestMappingsHistogramSeries);
        histogramOfMappingsByScore.createNewHistogramSeries(nameOfUniqueMappingsHistogramSeries);
        histogramOfMappingsByScore.createNewHistogramSeries(nameOfAllMappingsHistogramSeries);
        histogramOfMappingsByScore.createNewHistogramSeries(nameOfUniqueInclJunctionsHistogramSeries);

        TreeSet<ExtendedReadMapping> setOfExtendedReadMappingsSortedByPosition = new TreeSet<ExtendedReadMapping>(new ExtendedReadMappingByPositionComparator());

        //Counts end up in final report.
        int countOfReadsProcessed = 0;
        int countOfReadsWithAtLeastOneSplitHavingAFilterTag = 0;
        int countOfReadsWithAllSplitsHavingAFilterTag = 0;
        int countOfReadsWithTooFewMappings = 0;
        int countOfReadsWithTooManyMappings = 0;
        int countOfReadsWithNumberOfMappingsInRequiredRange = 0;

        int countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore = 0;
        int countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreThatAreUnique = 0;
        int countOfReadsWithReadSplitsAligningToTheSameLocation = 0;

        int countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[] = new int[filesMAXByRefPartitionByReadSequenceSplit[0].length];
        for (int i = 0; i < countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit.length; i++)
            countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[i] = 0;


        BufferedWriter writerMergedOutput = new BufferedWriter(new FileWriter(fileMergedOutput));
        BufferedWriter writerReadsWithoutMinScoringAlignment = new BufferedWriter(new FileWriter(fileReadsWithoutMinScoringAlignment));

        BufferedRandomAccessFile bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[][] = new BufferedRandomAccessFile[filesMAXByRefPartitionByReadSequenceSplit.length][filesMAXByRefPartitionByReadSequenceSplit[0].length];
        String lines[][] = new String[bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit.length][bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[0].length];
        
        JunctionFastaFile junctionFasta = null;
        File tmpDir = fileMergedOutput.getParentFile().getParentFile().getParentFile();
        if (fileJunctionMAX != null) {
        	junctionFasta = new JunctionFastaFile(new File(tmpDir, "junctions.fa"));
        }
        
        for (int indexRefFilePartition = 0; indexRefFilePartition < filesMAXByRefPartitionByReadSequenceSplit.length; indexRefFilePartition++) {
            for (int indexReadSequenceSplit = 0; indexReadSequenceSplit < filesMAXByRefPartitionByReadSequenceSplit[0].length; indexReadSequenceSplit++) {

                bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit] =
                                new BufferedRandomAccessFile(filesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit], "r");

                
                // advance to the (indexOfFirstReadToMerge +1)th read
                int indexOfCurrentRead = 0;
                File fileMAXFileIndex = new File(filesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].getPath() + extensionForMAXFileIndexFiles);
                if (fileMAXFileIndex.exists()) {
                	//Advance to the indexed read that appears before indexOfFirstReadToMerge
                    System.out.println("MAX file index file " + fileMAXFileIndex.getPath() + " found.");
                    TreeMap<String, Long> mapReadIndexToCharacterIndex = TextFileUtilities.loadStringToLongMapFromFile(fileMAXFileIndex, 0, 1, "\t", 0);

                    int indexOfReadToSeekTo = indexOfFirstReadToMerge;
                    long indexOfCharToSeekTo = -1;
                    while (indexOfCharToSeekTo < 0 && indexOfReadToSeekTo >= 0) {
                        if (mapReadIndexToCharacterIndex.containsKey("" + indexOfReadToSeekTo))
                            indexOfCharToSeekTo = mapReadIndexToCharacterIndex.get("" + indexOfReadToSeekTo);
                        else
                            indexOfReadToSeekTo--;
                    }    

                    bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].seek(indexOfCharToSeekTo);
//                    int numberOfCharsPerSkip = 2000000000;
//                    int numberOfSkips = (int)(indexOfCharToSeekTo / numberOfCharsPerSkip);
//                    for (int i = 0; i < numberOfSkips; i++) {
//                    	@SuppressWarnings("unused")
//                        int bytesSkipped = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].skipBytes(numberOfCharsPerSkip);
//                        //System.out.println("bytesSkipped:\t" + bytesSkipped);
//                    }
//                    @SuppressWarnings("unused")
//                    int bytesSkipped = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].skipBytes((int)(indexOfCharToSeekTo % numberOfCharsPerSkip));
//                    //System.out.println("bytesSkipped:\t" + bytesSkipped);
//                    //System.out.println("getFilePointer:\t" +  bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].getFilePointer());

                    
                    indexOfCurrentRead = indexOfReadToSeekTo;
                    lines[indexRefFilePartition][indexReadSequenceSplit] = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].readNextLine(); // header

                    System.out.println(indexRefFilePartition + " " + indexReadSequenceSplit + ":\tSkipped " + indexOfCurrentRead + " reads by seeking past " + indexOfCharToSeekTo + " characters.");
                    System.out.println(" Current read is:\t" + lines[indexRefFilePartition][indexReadSequenceSplit]);

                } else
                    System.out.println("MAX file index file " + fileMAXFileIndex.getPath() + " not found.  Advancing MAX fileReaders line-by-line... THIS MAY BE VERY SLOW!");
                
                //Read up until we hit indexOfFirstReadToMerge
                while (indexOfCurrentRead++ < indexOfFirstReadToMerge) {
                    bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].readNextLine();    // sequence
                    lines[indexRefFilePartition][indexReadSequenceSplit] = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].readNextLine(); // header
                    if (indexOfCurrentRead % 10000000 == 0)
                        System.out.println(indexRefFilePartition + " " + indexReadSequenceSplit + ":\tSkipped " + indexOfCurrentRead + " reads.");
                }

                System.out.println(" Current read is:\t" + lines[indexRefFilePartition][indexReadSequenceSplit]);
            }
        }
        
        BufferedRandomAccessFile brafJunctionMax = null;
        boolean mergingJunctionMappings = false;
        String junctionMaxLine = null;
        if (fileJunctionMAX != null) {
        	mergingJunctionMappings = true;
        	brafJunctionMax = new BufferedRandomAccessFile(fileJunctionMAX, "r");
        	File indexFile = new File(fileJunctionMAX.getParent(), fileJunctionMAX.getName()+WholeTranscriptomeAnalyzer.EXTENSION_FOR_INDEX_FILE);
        	int indexOfCurrentRead = 0;
        	if (indexFile != null) {
        		System.out.println("MAX file index file " + indexFile.getPath() + " found.");
                TreeMap<String, Long> mapReadIndexToCharacterIndex = TextFileUtilities.loadStringToLongMapFromFile(indexFile, 0, 1, "\t", 0);
                int indexOfReadToSeekTo = indexOfFirstReadToMerge;
                long indexOfCharToSeekTo = -1;
                while (indexOfCharToSeekTo < 0 && indexOfReadToSeekTo >= 0) {
                    if (mapReadIndexToCharacterIndex.containsKey("" + indexOfReadToSeekTo))
                        indexOfCharToSeekTo = mapReadIndexToCharacterIndex.get("" + indexOfReadToSeekTo);
                    else
                        indexOfReadToSeekTo--;
                }

                brafJunctionMax.seek(indexOfCharToSeekTo);

                indexOfCurrentRead = indexOfReadToSeekTo;
                junctionMaxLine = brafJunctionMax.readNextLine(); // header
                
                System.out.println("junction :\tSkipped " + indexOfCurrentRead + " reads by seeking past " + indexOfCharToSeekTo + " characters.");
                System.out.println(" Current read is:\t" + junctionMaxLine);

        	} else {
        		 System.out.println("MAX file index file " + fileJunctionMAX.getPath() + " not found.  Advancing MAX fileReaders line-by-line... THIS MAY BE VERY SLOW!");
        	}
            while (indexOfCurrentRead++ < indexOfFirstReadToMerge) {
                brafJunctionMax.readNextLine();    // sequence
                junctionMaxLine = brafJunctionMax.readNextLine(); // header
                if (indexOfCurrentRead % 10000000 == 0)
                    System.out.println("junctions :\tSkipped " + indexOfCurrentRead + " reads.");
            }
            System.out.println(" Current read is:\t" + junctionMaxLine);
        }
        
        String idOfRead = getBeadId(lines[0][0]);
        for (int i=0; i<lines.length; i++) {
        	for (int j=0; j<lines[i].length; j++) {
        		if (!getBeadId(lines[i][j]).equals(idOfRead))
        			throw new Exception("Max files are not aligned.");
        	}
        }
        if (mergingJunctionMappings && !getBeadId(junctionMaxLine).equals(idOfRead))
        	throw new Exception("Max files are not aligned.");
        System.out.println("Max files are aligned.  Proceeding with merge.");
       
        File partitionTableFile = new File(tmpDir, "reference.partition_table.out");
        Reader reader = null;
        Map<Integer, String> refNumToRefId;
        Map<String, Integer> refIdToRefNum;
        try {
        	reader = new BufferedReader(new FileReader(partitionTableFile));
        	refNumToRefId = getRefNumToRefId(reader);
        	reader.close();
        	reader = new BufferedReader(new FileReader(partitionTableFile));
        	refIdToRefNum = getRefIdToRefNum(reader);
        } finally {
        	if (reader != null) reader.close();
        }
        File fileMergedSortedGff = new File(fileMergedSortedOutput.getParentFile(), fileMergedSortedOutput.getName().replace(".csfasta",".gff"));
        PrintWriter gffWriter = new PrintWriter(new FileWriter(fileMergedSortedGff));
        
        
        // each loop iteration deals with one read
        while (lines[0][0] != null && countOfReadsProcessed < (indexOfLastReadToMerge - indexOfFirstReadToMerge +1)) {

            countOfReadsProcessed++;

            if (countOfReadsProcessed % 1000000 == 0)
                System.out.println("Processed " + countOfReadsProcessed + " reads.");

            int indexOfFirstComma = lines[0][0].indexOf(',');
            idOfRead = lines[0][0].substring(1);
            if (indexOfFirstComma > 0)
                idOfRead = lines[0][0].substring(1, indexOfFirstComma);
            idOfRead = getBeadId(idOfRead);

            boolean atLeastOneReadSplitHasAFilterTag = false;
            boolean allReadSplitsHaveFilterTag = true;
            boolean readSplitIsFiltered[] = new boolean[bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[0].length];
            boolean readSplitHasAlignmentWithMinAlignScore[] = new boolean[bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[0].length];
            boolean alignmentIsUnique = false;
            for (int i = 0; i < readSplitIsFiltered.length; i++) {
                readSplitIsFiltered[i] = lines[0][i].split(",")[0].endsWith(WholeTranscriptomeAnalyzer.SUFFIX_FOR_MARKING_READS_AS_FILTERED);
                allReadSplitsHaveFilterTag &= readSplitIsFiltered[i];
                atLeastOneReadSplitHasAFilterTag |= readSplitIsFiltered[i];
                readSplitHasAlignmentWithMinAlignScore[i] = false;
            }
            if (allReadSplitsHaveFilterTag)
                countOfReadsWithAllSplitsHavingAFilterTag++;
            if (atLeastOneReadSplitHasAFilterTag)
                countOfReadsWithAtLeastOneSplitHavingAFilterTag++;

            String sequenceOfRead = "";

            boolean foundMappingToSameLocationWithMinAlignScores = false;
            int countOfMappingsForThisRead = 0;
            ExtendedReadMappingSetForRead setOfSortedExtendedReadMappingsAll = new ExtendedReadMappingSetForRead(idOfRead, new ExtendedReadMappingByScoreComparator());
            ExtendedReadMappingSetForRead setOfSortedExtendedReadMappingsAboveMinAlignScore = new ExtendedReadMappingSetForRead(idOfRead, new ExtendedReadMappingByScoreComparator());
            for (int indexRefFilePartition = 0; indexRefFilePartition < bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit.length; indexRefFilePartition++) {

                ExtendedReadMapping readMappings[][] = new ExtendedReadMapping[bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition].length][];
                for (int indexReadSequenceSplit = 0; indexReadSequenceSplit < readMappings.length; indexReadSequenceSplit++) {

                    readMappings[indexReadSequenceSplit] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(lines[indexRefFilePartition][indexReadSequenceSplit].substring(1));

                    sequenceOfRead = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].readNextLine();
                    lines[indexRefFilePartition][indexReadSequenceSplit] = bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].readNextLine();
                }

                if (countOfMappingsForThisRead > maxMappingsAllowedBeforeNotReportingRead)
                    continue;   // no need to process all the alignments if the number of these already exceeds the limit (but do need to advance the readers as above)

                // compare all the read splits keeping only a single mapping when mappings correspond to the same location
                HashMap<String, ExtendedReadMapping> mapCollapsedMappingIdToMapping = new HashMap<String, ExtendedReadMapping>();
                for (int indexOfReadSequenceSplit = 0; indexOfReadSequenceSplit < readMappings.length; indexOfReadSequenceSplit++)
                    if (!readSplitIsFiltered[indexOfReadSequenceSplit]) {

                        for (int indexOfMapping = 0; indexOfMapping < readMappings[indexOfReadSequenceSplit].length; indexOfMapping++) {

                            if (readMappings[indexOfReadSequenceSplit][indexOfMapping].getScore() >= minAlignmentScoreForReportingAlignment)
                                readSplitHasAlignmentWithMinAlignScore[indexOfReadSequenceSplit] = true;

                            // FIXME: should convert the ref sequence index before creating the idCollapsedMapping with listsOfFullRefSequenceNumbersOrderedForEachReferencePartition

                            // alignments for split reads that map to same region of the genome will have the same idCollapsedMapping
                            String idCollapsedMapping = "";
                            if (readMappings[indexOfReadSequenceSplit][indexOfMapping].getPositionOfAlignmentStartInReferenceSequence() >= 0)
                                idCollapsedMapping = readMappings[indexOfReadSequenceSplit][indexOfMapping].getIndexOfMatchingReferenceSequence() + "_"
                                                    + (readMappings[indexOfReadSequenceSplit][indexOfMapping].getPositionOfAlignmentStartInReferenceSequence() - readMappings[indexOfReadSequenceSplit][indexOfMapping].getPositionOfAlignmentStartInRead() +1);
                            else
                                idCollapsedMapping = readMappings[indexOfReadSequenceSplit][indexOfMapping].getIndexOfMatchingReferenceSequence() + "_"
                                                    + (readMappings[indexOfReadSequenceSplit][indexOfMapping].getPositionOfAlignmentStartInReferenceSequence()
                                                        + (sequenceOfRead.length() - readMappings[indexOfReadSequenceSplit][indexOfMapping].getLengthOfAlignment() - readMappings[indexOfReadSequenceSplit][indexOfMapping].getPositionOfAlignmentStartInRead()));

                            if (mapCollapsedMappingIdToMapping.containsKey(idCollapsedMapping)) { // keep the higher scoring alignment     FIXME: is this what we want?

                                int scoreOfExistingMapping = mapCollapsedMappingIdToMapping.get(idCollapsedMapping).getScore();
                                if (readMappings[indexOfReadSequenceSplit][indexOfMapping].getScore() > scoreOfExistingMapping)
                                    mapCollapsedMappingIdToMapping.put(idCollapsedMapping, readMappings[indexOfReadSequenceSplit][indexOfMapping]);

                                if (scoreOfExistingMapping >= minAlignmentScoreForReportingAlignment &&
                                        readMappings[indexOfReadSequenceSplit][indexOfMapping].getScore() >= minAlignmentScoreForReportingAlignment)
                                    foundMappingToSameLocationWithMinAlignScores = true;
                            } else
                                mapCollapsedMappingIdToMapping.put(idCollapsedMapping, readMappings[indexOfReadSequenceSplit][indexOfMapping]);
                        }
                    }

                Iterator<String> iteratorCollapsedMappingIds = mapCollapsedMappingIdToMapping.keySet().iterator();
                while (iteratorCollapsedMappingIds.hasNext()) {
                    String idCollapsedMapping = iteratorCollapsedMappingIds.next();
                    ExtendedReadMapping extendedReadMapping = mapCollapsedMappingIdToMapping.get(idCollapsedMapping);
                    countOfMappingsForThisRead++;

                    if (listsOfFullRefSequenceNumbersOrderedForEachReferencePartition != null) {    // re-assign the reference sequence index in accordance with the pre-split file
                    	extendedReadMapping.setIndexOfMatchingReferenceSequence(
                                listsOfFullRefSequenceNumbersOrderedForEachReferencePartition[indexRefFilePartition].get(extendedReadMapping.getIndexOfMatchingReferenceSequence() -1));    // ?????
                        extendedReadMapping.setIdOfMatchingReferenceSequence(extendedReadMapping.getIndexOfMatchingReferenceSequence() + "");
                    }

                    setOfSortedExtendedReadMappingsAll.addMapping(extendedReadMapping);

                    if (extendedReadMapping.getScore() >= minAlignmentScoreForReportingAlignment)
                        setOfSortedExtendedReadMappingsAboveMinAlignScore.addMapping(extendedReadMapping);

                }
            }


            if (!allReadSplitsHaveFilterTag && (!atLeastOneReadSplitHasAFilterTag || !filterReadsWithOneOrMoreSplitsHavingAFilterTag)) {

                if (countOfMappingsForThisRead >= minMappingsRequiredBeforeReportingRead && countOfMappingsForThisRead <= maxMappingsAllowedBeforeNotReportingRead) {

                    countOfReadsWithNumberOfMappingsInRequiredRange++;

                    Iterator<ExtendedReadMapping> iteratorOverExtendedReadMappings = setOfSortedExtendedReadMappingsAll.getSetOfSortedReadMappings().iterator();
                    ExtendedReadMapping extendedReadMappingBestFound = null;
                    ExtendedReadMapping extendedReadMappingSecondBestFound = null;
                    if (iteratorOverExtendedReadMappings.hasNext())
                        extendedReadMappingBestFound = iteratorOverExtendedReadMappings.next();
                    if (iteratorOverExtendedReadMappings.hasNext())
                        extendedReadMappingSecondBestFound = iteratorOverExtendedReadMappings.next();

                    
                    if ((extendedReadMappingSecondBestFound == null && (extendedReadMappingBestFound.getScore() - minScoreGapToSecondBestAlignmentForUniqueness) >= minPossibleAlignmentScore -1)
                            || (extendedReadMappingSecondBestFound != null && (extendedReadMappingBestFound.getScore() - extendedReadMappingSecondBestFound.getScore()) >= minScoreGapToSecondBestAlignmentForUniqueness)) {
                        histogramOfMappingsByScore.addObservationToHistogramSeries(nameOfUniqueMappingsHistogramSeries, extendedReadMappingBestFound.getScore());
                        alignmentIsUnique = true;
                    }
                    if (extendedReadMappingBestFound != null)
                        histogramOfMappingsByScore.addObservationToHistogramSeries(nameOfBestMappingsHistogramSeries, extendedReadMappingBestFound.getScore());

                    Iterator<ExtendedReadMapping> extendedReadMappingsAllIterator = setOfSortedExtendedReadMappingsAll.getSetOfSortedReadMappings().iterator();
                    while (extendedReadMappingsAllIterator.hasNext())
                        histogramOfMappingsByScore.addObservationToHistogramSeries(nameOfAllMappingsHistogramSeries, extendedReadMappingsAllIterator.next().getScore());

                    if (foundMappingToSameLocationWithMinAlignScores) countOfReadsWithReadSplitsAligningToTheSameLocation++;

                    if (setOfSortedExtendedReadMappingsAboveMinAlignScore.size() > 0) {

                        countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore++;

                        for (int i = 0; i < readSplitHasAlignmentWithMinAlignScore.length; i++) {
                            if (readSplitHasAlignmentWithMinAlignScore[i])
                                countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[i]++;
                        }

                        writerMergedOutput.write(">" + setOfSortedExtendedReadMappingsAboveMinAlignScore.toString());
                        writerMergedOutput.newLine();
                        writerMergedOutput.write(sequenceOfRead);
                        writerMergedOutput.newLine();

                        if (alignmentIsUnique) {

                            countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreThatAreUnique++;

                            extendedReadMappingBestFound.setSequenceOfRead(sequenceOfRead);
                            setOfExtendedReadMappingsSortedByPosition.add(extendedReadMappingBestFound);

                        }

                    }

                    //Combine Junction Mappings with reads having min score.
    	            List<ExtendedReadMapping> allMappings = new ArrayList<ExtendedReadMapping>(setOfSortedExtendedReadMappingsAll.getSetOfSortedReadMappings());
    	            if (brafJunctionMax != null) {
    	            	
    	            	//Combine merged genomic max with junction max {
    	            	ExtendedReadMapping[] junctionMappings = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(junctionMaxLine.substring(1));
    	           			
    	        		for (ExtendedReadMapping mapping : junctionMappings) {
    	        			if (mapping.getPositionOfAlignmentStartInReferenceSequence() < 0) continue; //Only positive strand alignments.
    	        			FlankedJunction junction = junctionFasta.getFlankedJunction(mapping.getIndexOfMatchingReferenceSequence());
    	        			if (junction.isKnown() ) 
    	        				mapping.setScore(mapping.getScore() - penaltyForKnownJunctions);
    	        			else
    	        				mapping.setScore(mapping.getScore() - penaltyForPutativeJunctions);
    	            		JunctionMapping jMapping= new JunctionMapping(mapping, junction, refIdToRefNum.get(junction.getChromId()));
    	            		//Only add mapping if it actually spans the junction and satisfies the minimum score.
    	            		if (jMapping.isJunctionSpanned() &&
    	            			jMapping.getScore() > minAlignmentScoreForReportingAlignment)
    	            			allMappings.add(jMapping);        		
    	            	}
    	            } 
    	            
    	    		//Sort the mappings by score
    				Collections.sort(allMappings, ExtendedReadMappingByScoreComparator.INSTANCE);
    				
    				//Collapse until first two mappings are the junction and genomic mappings to the same location.
    				while (true) {
    					if (allMappings.size() < 2) break;
    					int indexOfJMapping = -1;
    					int indexOfEMapping = -1;
    					if (allMappings.get(0) instanceof JunctionMapping) indexOfJMapping = 0;
    					else indexOfEMapping = 0;
    					if (allMappings.get(1) instanceof JunctionMapping) indexOfJMapping = 1;
    					else indexOfEMapping = 1;
    	
    					if (indexOfJMapping > -1 && indexOfEMapping > -1) {
    						JunctionMapping jMapping = (JunctionMapping)allMappings.get(indexOfJMapping);
    						ExtendedReadMapping eMapping = allMappings.get(indexOfEMapping);
    						if (jMapping.mapsToSameLocationAs(eMapping)) {
    							//Select the mapping with the higher score.
    							if (jMapping.getScore() > eMapping.getScore()) {
    								allMappings.remove(indexOfEMapping);
    							} else {
    								allMappings.remove(indexOfJMapping);
    							}
    						} else {
    							//First two mappings do not map to same location.
    							break;
    						}
    					} else {
    						//First two mappings are same type.
    						break;
    					}
    				}
    				
    				//Determine the unique mapping.
    				ExtendedReadMapping uniqueMapping = null;
    				if (allMappings.size() < minMappingsRequiredBeforeReportingRead) {
    					uniqueMapping = null;
    				} else if (allMappings.size() == 1) {
    					if ((allMappings.get(0).getScore() - minScoreGapToSecondBestAlignmentForUniqueness) >= minPossibleAlignmentScore -1 )
    						uniqueMapping = allMappings.get(0);
    				} else if (allMappings.size() <= maxMappingsAllowedBeforeNotReportingRead &&
    						 allMappings.get(0).getScore() - allMappings.get(1).getScore() >= minScoreGapToSecondBestAlignmentForUniqueness) {
    					uniqueMapping = allMappings.get(0);
    				}
    	
    				if (alignmentIsUnique != (uniqueMapping != null)) {
    					System.out.printf("Uniqueness discrepancy: %s %b %s %s %s\n", idOfRead, alignmentIsUnique, allMappings, countOfMappingsForThisRead, maxMappingsAllowedBeforeNotReportingRead);
    				}
    	
    	    		//If it is unique, create GFF features.
    	    		List<GffFeature> features = new ArrayList<GffFeature>();
    	    		String gffAttributes = "bd=%s;rs=%d;mm=%d;g=%s;i=%s;";
    	    		if (uniqueMapping != null) {
    	    			if (uniqueMapping instanceof JunctionMapping) {
    	    				//Junction mappings are represented as two entries in the GFF file,
    	    				//representing the portions of reads that map to each size of the junction.
    	    				JunctionMapping jMapping = (JunctionMapping)uniqueMapping;
    	    				GffFeature leftFeature =
    	    					new GffFeature(jMapping.getJunction().getChromId(),
    	    							       jMapping.getStrand(),
    	    								   jMapping.getLeftStart(),
    	    								   jMapping.getLeftEnd());
    	    				String attributes = String.format(gffAttributes, idOfRead, uniqueMapping.getPositionOfAlignmentStartInRead(), uniqueMapping.getNumberOfMismatches(), sequenceOfRead, refIdToRefNum.get(jMapping.getJunction().getChromId()));
    	    				leftFeature.setSource("wtp");
    	    				leftFeature.setAdditionalInfo(attributes + String.format("jp=%d;", jMapping.getRightStart()) + String.format("jt=%s;", jMapping.getJunction().isKnown() ? "k" : "p"));
    	    				leftFeature.setFeatureType("read");
    	    				leftFeature.setScore(uniqueMapping.getScore());
    	    				features.add(leftFeature);
    	    				GffFeature rightFeature =
    	    					new GffFeature(jMapping.getJunction().getChromId(),
    							       jMapping.getStrand(),
    								   jMapping.getRightStart(),
    								   jMapping.getRightEnd());
    	    				rightFeature.setAdditionalInfo(attributes + String.format("jp=%d;", jMapping.getLeftStart()) + String.format("jt=%s;", jMapping.getJunction().isKnown() ? "k" : "p"));
    	    				rightFeature.setSource("wtp");
    	    				rightFeature.setFeatureType("read");
    	    				rightFeature.setScore(uniqueMapping.getScore());
    	    				features.add(rightFeature);
    	    			} else {
    	    				//Genomic mappings are represented as single entries in the GFF file.
    	    				String seqId = refNumToRefId.get(uniqueMapping.getIndexOfMatchingReferenceSequence()).toString();
    	    				GffFeature feature = 
    	    					new GffFeature(seqId, 
    	    							uniqueMapping.getPositionOfAlignmentStartInReferenceSequence() < 0 ? Strand.NEGATIVE : Strand.POSITIVE,
    	    							Math.abs(uniqueMapping.getPositionOfAlignmentStartInReferenceSequence())+1,
    	    							Math.abs(uniqueMapping.getPositionOfAlignmentStartInReferenceSequence())+uniqueMapping.getLengthOfAlignment());
    	    				feature.setAdditionalInfo(String.format(gffAttributes, idOfRead, uniqueMapping.getPositionOfAlignmentStartInRead(), uniqueMapping.getNumberOfMismatches(), sequenceOfRead, uniqueMapping.getIndexOfMatchingReferenceSequence()));
    	    				feature.setSource("wtp");
    	    				feature.setFeatureType("read");
    	    				feature.setScore(uniqueMapping.getScore());
    	    				features.add(feature);
    	    			}
    	    			
    	    			histogramOfMappingsByScore.addObservationToHistogramSeries(nameOfUniqueInclJunctionsHistogramSeries, uniqueMapping.getScore());
    	    		}
    	    		
    	    		//Write the gff features to the file.
    	    		for (GffFeature feature : features) {
    	    			feature.toGffRow(gffWriter);
    	    		}
    	    		if (brafJunctionMax != null) {
    	    			brafJunctionMax.readNextLine(); //Discard Sequence line
    	    			junctionMaxLine = brafJunctionMax.readNextLine(); //Header for next loop.
    	    		}
                    
                } else { //end if (countOfMappingsForThisRead >= minMappingsRequiredBeforeReportingRead && countOfMappingsForThisRead <= maxMappingsAllowedBeforeNotReportingRead)

                    if (countOfMappingsForThisRead < minMappingsRequiredBeforeReportingRead)
                        countOfReadsWithTooFewMappings++;

                    if (countOfMappingsForThisRead > maxMappingsAllowedBeforeNotReportingRead)
                        countOfReadsWithTooManyMappings++;
                }

                if (setOfSortedExtendedReadMappingsAboveMinAlignScore.size() == 0) {
                    writerReadsWithoutMinScoringAlignment.write(">" + idOfRead);
                    writerReadsWithoutMinScoringAlignment.newLine();
                    writerReadsWithoutMinScoringAlignment.write(sequenceOfRead);
                    writerReadsWithoutMinScoringAlignment.newLine();
                }
            } // End if (!allReadSplitsHaveFilterTag && (!atLeastOneReadSplitHasAFilterTag || !filterReadsWithOneOrMoreSplitsHavingAFilterTag))
        } // Done loop per bead ID.
        
        for (int indexRefFilePartition = 0; indexRefFilePartition < filesMAXByRefPartitionByReadSequenceSplit.length; indexRefFilePartition++)
            for (int indexReadSequenceSplit = 0; indexReadSequenceSplit < filesMAXByRefPartitionByReadSequenceSplit[0].length; indexReadSequenceSplit++)
                bufferedRandomAccessFilesMAXByRefPartitionByReadSequenceSplit[indexRefFilePartition][indexReadSequenceSplit].close();

        writerMergedOutput.close();
        writerReadsWithoutMinScoringAlignment.close();
        if (junctionFasta != null)  junctionFasta.close();
        if (gffWriter != null) gffWriter.close();

        BufferedWriter writerMergedSortedOutput = new BufferedWriter(new FileWriter(fileMergedSortedOutput));
        Iterator<ExtendedReadMapping> iteratorOverSortedReadMappings = setOfExtendedReadMappingsSortedByPosition.iterator();
        while (iteratorOverSortedReadMappings.hasNext()) {
            ExtendedReadMapping extendedReadMapping = iteratorOverSortedReadMappings.next();
            writerMergedSortedOutput.write(">" + extendedReadMapping.getIdOfRead() + extendedReadMapping.toString());
            writerMergedSortedOutput.newLine();
            writerMergedSortedOutput.write(extendedReadMapping.getSequenceOfRead());
            writerMergedSortedOutput.newLine();
        }
        writerMergedSortedOutput.close();
        
        //Sort the gffFile
        FileReader readerMergedSortedGff = null;
        try {
        	readerMergedSortedGff = new FileReader(fileMergedSortedGff);
        	LineUtils.sortLines(readerMergedSortedGff, fileMergedSortedGff, GffFeature.FORMAT);
        } finally {
        	if (readerMergedSortedGff != null) readerMergedSortedGff.close();
        }
        
        reportOfExtendedReadMappings.setMinMappingsRequiredBeforeReportingRead(minMappingsRequiredBeforeReportingRead);
        reportOfExtendedReadMappings.setMaxMappingsAllowedBeforeNotReportingRead(maxMappingsAllowedBeforeNotReportingRead);
        reportOfExtendedReadMappings.setMinAlignmentScoreForReportingMapping(minAlignmentScoreForReportingAlignment);
        reportOfExtendedReadMappings.setMinScoreGapToSecondBestAlignmentForUniqueness(minScoreGapToSecondBestAlignmentForUniqueness);
        reportOfExtendedReadMappings.setCountOfReadsProcessed(countOfReadsProcessed);
        reportOfExtendedReadMappings.setCountOfReadsWithAtLeastOneSplitHavingAFilterTag(countOfReadsWithAtLeastOneSplitHavingAFilterTag);
        reportOfExtendedReadMappings.setCountOfReadsWithAllSplitsHavingAFilterTag(countOfReadsWithAllSplitsHavingAFilterTag);
        reportOfExtendedReadMappings.setCountOfReadsWithTooFewMappings(countOfReadsWithTooFewMappings);
        reportOfExtendedReadMappings.setCountOfReadsWithTooManyMappings(countOfReadsWithTooManyMappings);
        reportOfExtendedReadMappings.setCountOfReadsWithNumberOfMappingsInRequiredRange(countOfReadsWithNumberOfMappingsInRequiredRange);
        reportOfExtendedReadMappings.setCountOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore(countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore);
        reportOfExtendedReadMappings.setCountOfReadsUniquelyAlignedWithMinAlignScore(countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreThatAreUnique);
        reportOfExtendedReadMappings.setCountsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit(countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit);
        reportOfExtendedReadMappings.setCountOfReadsWithReadSplitsAligningToTheSameLocation(countOfReadsWithReadSplitsAligningToTheSameLocation);
        reportOfExtendedReadMappings.setHistogramOfMappingsByScore(histogramOfMappingsByScore);


    }

    private static String getBeadId(String s) {
    	if (s == null ) return null;
    	s = s.replaceAll(",.*", ""); //strip hits.
    	return s.replaceAll(WholeTranscriptomeAnalyzer.SUFFIX_FOR_MARKING_READS_AS_FILTERED + "$", ""); //Strip Filtered label.
    }
    
    private static Map<Integer, String> getRefNumToRefId(Reader reader) throws IOException {
    	BufferedReader br = new BufferedReader(reader);
    	Map<Integer, String> map = new HashMap<Integer, String>();
    	for (String line = br.readLine(); line != null; line = br.readLine()) {
    		if (line.trim().startsWith("#")) continue;
    		String[] fields = line.split("\t");
    		Integer sequenceNumberInRefPartition = Integer.valueOf(fields[2]);
    		String sequenceId = fields[3].replaceAll("\\s.*$", "");
    		map.put(sequenceNumberInRefPartition, sequenceId);
    		//map.put(sequenceId, sequenceNumberInRefPartition);
    	}
    	return map;
    }
    
    private static Map<String, Integer> getRefIdToRefNum(Reader reader) throws IOException {
    	BufferedReader br = new BufferedReader(reader);
    	Map<String, Integer> map = new HashMap<String, Integer>();
    	for (String line = br.readLine(); line != null; line = br.readLine()) {
    		if (line.trim().startsWith("#")) continue;
    		String[] fields = line.split("\t");
    		Integer sequenceNumberInRefPartition = Integer.valueOf(fields[2]);
    		String sequenceId = fields[3].replaceAll("\\s.*$", "");
    		map.put(sequenceId, sequenceNumberInRefPartition);
    	}
    	return map;
    }
}

class JunctionMapping extends ExtendedReadMapping {

	private FlankedJunction junction;
    private int minNumberOfColorsInAFlank = 4;
    private int indexOfMatchingReferenceSequence;
    
	public JunctionMapping(ExtendedReadMapping e, FlankedJunction junction, int indexOfReferenceSequence) {
		super(e.getIdOfRead(), junction.getChromId(), e.getIndexOfMatchingReferenceSequence(), e.getPositionOfAlignmentStartInReferenceSequence(), e.getNumberOfMismatches(), e.getScore(), e.getPositionOfAlignmentStartInRead(), e.getLengthOfAlignment());
		this.junction = junction;
		this.indexOfMatchingReferenceSequence = indexOfReferenceSequence;
	}
	
	@Override
	public String toString() {
		return junction.toString() + " " + super.toString();
	}
	
	@Override
	public int getIndexOfMatchingReferenceSequence() {
		return indexOfMatchingReferenceSequence;
	}
	
	public FlankedJunction getJunction() {
		return this.junction;
	}

	public Strand getStrand() {
		return this.junction.getStrand();
	}
	
	public int getNumColorsInDonorSequence() {
		if (this.getPositionOfAlignmentStartInReferenceSequence() > junction.getLastPositionInDonorSequence()) {
			return 0;
		}
		return junction.getLastPositionInDonorSequence() - this.getPositionOfAlignmentStartInReferenceSequence() + 1;
	}
	
	public int getNumColorsInAcceptorSequence() {
		return this.getLengthOfAlignment() - getNumColorsInDonorSequence();
	}
	
	public boolean isJunctionSpanned() {
		//
		if (getNumColorsInDonorSequence() < minNumberOfColorsInAFlank ||
			getNumColorsInAcceptorSequence() < minNumberOfColorsInAFlank) return false;
		return true;
	}
	
	public int getLeftStart() {
		if (this.junction.getStrand() == Strand.POSITIVE ) {
			return this.junction.getPositionInGenomicReference(this.getPositionOfAlignmentStartInReferenceSequence()) + 1;
		} else {
			return this.junction.getPositionInGenomicReference(this.getPositionOfAlignmentStartInReferenceSequence() + this.getLengthOfAlignment() - 1) + 1;
		}
	}
	
	public int getLeftEnd() {
		return this.junction.getLeftEnd();
	}
	
	public int getRightStart() {
		return this.junction.getRightStart();
	}
	
	public int getRightEnd() {
		if (this.junction.getStrand() == Strand.POSITIVE) {
			return this.junction.getPositionInGenomicReference(this.getPositionOfAlignmentStartInReferenceSequence() + this.getLengthOfAlignment() - 1) + 1;
		} else {
			return this.junction.getPositionInGenomicReference(this.getPositionOfAlignmentStartInReferenceSequence()) + 1;
		}
	}
	
	//Return true if the extended read mapping maps to one side of the junction or the other.
	public boolean mapsToSameLocationAs(ExtendedReadMapping that) {
		if (!that.getClass().equals(ExtendedReadMapping.class)) throw new IllegalArgumentException("Argument must be class ExtendedReadMapping, subclasses not supported.");
		Strand thatStrand = that.getPositionOfAlignmentStartInReferenceSequence() >=0 ? Strand.POSITIVE : Strand.NEGATIVE;
		if (this.getIndexOfMatchingReferenceSequence() != that.getIndexOfMatchingReferenceSequence()) return false;
		if (this.getStrand() != thatStrand) return false;
		int thatLeftmostPosition = Math.abs(that.getPositionOfAlignmentStartInReferenceSequence());
		if (thatLeftmostPosition >= this.junction.getLeftStart() - 1 && thatLeftmostPosition <= this.junction.getLeftEnd() - 1) return true;
		int thatRightmostPosition = Math.abs(that.getPositionOfAlignmentStartInReferenceSequence()) + that.getLengthOfAlignment() - 1;
		if (thatRightmostPosition >= this.junction.getRightStart() - 1 && thatRightmostPosition <= junction.getRightEnd() - 1) return true;
		return false;
	}
	
	public static void main(String[] args) {
		FlankedJunction junction = FlankedJunction.parseFlankedJunction("junction:1:2073493-2073538|2073911-2073956:foo:+:known");
	    JunctionMapping jMapping = new JunctionMapping(new ExtendedReadMapping("read_1", "1", 1, 43, 1, 42, 1, 45), junction, 1);
	    System.out.println(jMapping.getNumColorsInDonorSequence());
	    System.out.println(jMapping.getNumColorsInAcceptorSequence());
	    System.out.println(jMapping.isJunctionSpanned());
	    System.out.printf("%d-%d\n",jMapping.getLeftStart(), jMapping.getLeftEnd());
	    System.out.printf("%d-%d\n",jMapping.getRightStart(), jMapping.getRightEnd());
	}
}

