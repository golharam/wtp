package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.mapper.FilteringMode;
import com.lifetechnologies.solid.wt.mapper.FilteringTask;


/**
 * User: tuchbb
 * Date: Sep 26, 2008
 * Time: 12:30:54 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class WholeTranscriptomeAnalyzer {

	public static final String WT_HOME_SYSTEM_PROPERTY = "com.lifetechnologies.solid.wt.home";
	public static final String EXPERIMENTAL_MODE_SYSTEM_PROPERTY = "com.lifetechnologies.solid.wt.experimentalMode";

	private static PropertyResourceBundle resourceBundle;
	
    private static String VERSION;
    private static String NAME_OF_TOOL;

    private static String SCHEDULING_ENVIRONMENT = "pbs";
    private static String NAME_OF_QUEUE;
    private static Long MAX_MEMORY_PER_JOB_IN_BYTES;
    private static Long MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES;
    private static Integer NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE;
    private static Integer NUMBER_OF_READS_FILE_PARTITIONS = 1; // FIXME: do not change this value until the merge of readsFilePartitions is updated to also merge index files
    private static Boolean JOBS_RERUNNABLE;
    private static Float MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR = 1.0f;
    private static String SCHEDULER_RESOURCE_REQUIREMENTS;

    private static String ADDITIONAL_SCHEDULER_OPTIONS;
    
    private static Boolean RUN_READ_SPLITTING;
    private static Boolean RUN_READ_FILTERING;
    private static Boolean RUN_REFERENCE_PARTITIONING;
    private static Boolean RUN_MAPPING;
    private static Boolean RUN_EXTENSION;
    private static Boolean RUN_MERGE;
    

    private static FilteringMode FILTERING_MODE;
    
    private static Boolean COMPRESS_INTERMEDIATE_FILES;
    private static Boolean DELETE_INTERMEDIATE_FILES;

    private static Integer LENGTH_OF_READS;
    private static String MASK_OF_READS;
    private static Integer MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ;
    private static Integer MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ;
    private static Boolean VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION;
    private static Boolean MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION;
    private static Integer LENGTH_OF_FIRST_PART_OF_READ;
    private static Integer LENGTH_OF_LAST_PART_OF_READ;
    private static Integer MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING;
    private static Integer MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING;
    private static Integer MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING;
    private static Integer MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT;
    private static Integer MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS;
    private static Integer SCORE_OF_MATCH_FOR_EXTENSION;
    private static Integer SCORE_OF_MISMATCH_FOR_EXTENSION;

    private static File FOLDER_INSTALLATION_ROOT;
    private static File FILE_FILTER_SEQUENCES_FASTA;
    private static File FILE_REFERENCE_FASTA;
    private static File FILE_FULL_LENGTH_READS_CSFASTA;
    private static File FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES;
    private static File FOLDER_FOR_OUTPUT_FILES;

    public static String EXTENSION_FOR_READS_FILE = ".csfasta";
    public static String EXTENSION_FOR_MAPREADS_OUTPUT_FILE = ".ma";
    public static String EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE = ".max";
    public static String EXTENSION_FOR_MERGED_OUTPUT_FILE = ".merged";
    public static String EXTENSION_FOR_FILTERED_OUTPUT_FILE = ".filtered";
    public static String EXTENSION_FOR_SORTED_OUTPUT_FILE = ".sorted";
    public static String EXTENSION_FOR_READS_WITHOUT_MIN_SCORING_ALIGNMENT_FILE = ".readsWithoutMinScoringAlignment";
    public static String EXTENSION_FOR_COMPRESSED_FILE = ".zip";
    public static String EXTENSION_FOR_INDEX_FILE = ".index";
    public static String SUFFIX_FOR_MARKING_READS_AS_FILTERED = "_FILT";

    public static File FOLDER_CONTAINING_SCHEMA_FILES;
    public static File WHOLE_TRANSCRIPTOME_JAR;
    
    public static final FilenameFilter FILENAME_FILTER_MAX = new FilenamePrefixSuffixFilter("", EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
    public static final FilenameFilter FILENAME_FILTER_MERGED_MAX = new FilenamePrefixSuffixFilter("", EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
    public static final FilenameFilter FILENAME_FILTER_MERGED_MAX_COMPRESSED = new FilenamePrefixSuffixFilter("", EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + EXTENSION_FOR_COMPRESSED_FILE);

    /**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		loadResourceBundle();
		NAME_OF_TOOL = resourceBundle.getString("application.name");
		VERSION = resourceBundle.getString("version");
       
        
        if (args.length < 1 || args[0].equalsIgnoreCase("-v")) {
            System.out.println(NAME_OF_TOOL + " " + VERSION);
            return;
        }
        
        System.out.println("Starting "+NAME_OF_TOOL+", Version: "+VERSION+" at " + new Date());
        
        //Get Installation dir from the system property set on command line with -D
        FOLDER_INSTALLATION_ROOT = new File(System.getProperty(WT_HOME_SYSTEM_PROPERTY));
        if (FOLDER_INSTALLATION_ROOT.exists() == false)
        	throw new Exception(FOLDER_INSTALLATION_ROOT + "Does not exist.");
        WHOLE_TRANSCRIPTOME_JAR = new File(FOLDER_INSTALLATION_ROOT, "pkg/WholeTranscriptome.jar");
        if (WHOLE_TRANSCRIPTOME_JAR.exists() == false)
        	throw new Exception(FOLDER_INSTALLATION_ROOT + "is not a Whole Transcriptome Pipeline folder");
        
        FOLDER_CONTAINING_SCHEMA_FILES = new File(FOLDER_INSTALLATION_ROOT.getPath() + "/etc/schemas/");
        if (!FOLDER_CONTAINING_SCHEMA_FILES.exists())
            throw new Exception("Could not find folder containing schema files at: " + FOLDER_CONTAINING_SCHEMA_FILES.getPath());
        
        loadConfigurationFromFile(new File(args[0]));
        
        // this gets around a bug associated with the GUI invoked by JFreeChart when X is no available
        System.setProperty("java.awt.headless", "true"); 

        File fileMapReadsExe = new File(FOLDER_INSTALLATION_ROOT.getPath() + "/pkg/mapreads");   // please ensure that this is version 1.1.6 or later
        if (!fileMapReadsExe.exists())
            throw new Exception("Could not find mapReads exe at: " + fileMapReadsExe.getPath());

        
        File fileExtendMapReadsExe = new File(FOLDER_INSTALLATION_ROOT.getPath() + "/pkg/extendMappedReads.py");
        if (!fileExtendMapReadsExe.exists())
            throw new Exception("Could not find extendMappedReads exe at: " + fileExtendMapReadsExe.getPath());
        
        //File folderContainingJavaClasses = new File(FOLDER_INSTALLATION_ROOT.getPath() + "/bin/");
        //if (!folderContainingJavaClasses.exists())
        //    throw new Exception("Could not find folder containing java classes at: " + folderContainingJavaClasses.getPath());

        File folderContainingJavaLibraries = new File(FOLDER_INSTALLATION_ROOT.getPath() + "/lib/");
            if (!folderContainingJavaLibraries.exists())
                throw new Exception("Could not find folder containing java libraries at: " + folderContainingJavaLibraries.getPath());


        File folderForTempFiles = new File(FOLDER_FOR_OUTPUT_FILES.getPath() + "/tmp/");
        File folderForFilterMapping = new File(folderForTempFiles.getPath() + "/filtering/");
        File folderForMappingAndExtension = new File(folderForTempFiles.getPath() + "/mappingAndExtension/");
        File folderForMerge = new File(folderForTempFiles.getPath() + "/merge/");
        File folderForScriptFiles = new File(FOLDER_FOR_OUTPUT_FILES.getPath() + "/scripts/");
        File folderForImportantOutputFiles = new File(FOLDER_FOR_OUTPUT_FILES.getPath() + "/output/"); // Important = those likely to be kept for posterity
        File fileReferencePartitionTable = new File(folderForTempFiles.getPath() + "/reference.partition_table.out");
        File fileReadsCount = new File(folderForTempFiles.getPath() + "/reads.count.out");

        if (!FOLDER_FOR_OUTPUT_FILES.exists())    FOLDER_FOR_OUTPUT_FILES.mkdir();
        if (!folderForTempFiles.exists())    folderForTempFiles.mkdir();
        if (!folderForScriptFiles.exists())  folderForScriptFiles.mkdirs();
        if (!folderForImportantOutputFiles.exists())  folderForImportantOutputFiles.mkdirs();

        File fileMasterFilteringJobSubmissionScript = new File(folderForScriptFiles.getPath() + "/filtering.master.script_submission.sh");
        File fileMasterFilteringJobRemovalScript = new File(folderForScriptFiles.getPath() + "/filtering.master.script_removal.sh");
        File fileMasterFilteringListOfScriptOutputFiles = new File(folderForTempFiles.getPath() + "/filtering.list_of_script_output_files.out");

        File fileMasterMappingJobSubmissionScript = new File(folderForScriptFiles.getPath() + "/mapping.master.script_submission.sh");
        File fileMasterMappingJobRemovalScript = new File(folderForScriptFiles.getPath() + "/mapping.master.script_removal.sh");
        File fileMasterMappingListOfScriptOutputFiles = new File(folderForTempFiles.getPath() + "/mapping.list_of_script_output_files.out");

        File fileMasterExtensionJobSubmissionScript = new File(folderForScriptFiles.getPath() + "/extension.master.script_submission.sh");
        File fileMasterExtensionJobRemovalScript = new File(folderForScriptFiles.getPath() + "/extension.master.script_removal.sh");
        File fileMasterExtensionListOfScriptOutputFiles = new File(folderForTempFiles.getPath() + "/extension.list_of_script_output_files.out");

        File fileMasterMergeJobSubmissionScript = new File(folderForScriptFiles.getPath() + "/merge.master.script_submission.sh");
        File fileMasterMergeJobRemovalScript = new File(folderForScriptFiles.getPath() + "/merge.master.script_removal.sh");
        File fileMasterMergeListOfScriptOutputFiles = new File(folderForTempFiles.getPath() + "/merge.list_of_script_output_files.out");


        File fileFullyMergedMAXcsfasta = new File(folderForImportantOutputFiles.getPath() + "/" + FILE_FULL_LENGTH_READS_CSFASTA.getName().substring(0, FILE_FULL_LENGTH_READS_CSFASTA.getName().lastIndexOf('.'))
                                                    + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE);
        File fileFullyMergedAndSortedMAXcsfasta = new File(folderForImportantOutputFiles.getPath() + "/" + FILE_FULL_LENGTH_READS_CSFASTA.getName().substring(0, FILE_FULL_LENGTH_READS_CSFASTA.getName().lastIndexOf('.'))
                                                    + EXTENSION_FOR_SORTED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE);
        File fileReadsWithoutMinScoringAlignmentCSfasta = new File(folderForImportantOutputFiles.getPath() + "/" + FILE_FULL_LENGTH_READS_CSFASTA.getName().substring(0, FILE_FULL_LENGTH_READS_CSFASTA.getName().lastIndexOf('.'))
                                                    + EXTENSION_FOR_READS_WITHOUT_MIN_SCORING_ALIGNMENT_FILE + EXTENSION_FOR_READS_FILE);

        File fileFilterReportTxt = new File(folderForImportantOutputFiles.getPath() + "/filterReport.txt");
        File fileMappingReportTxt = new File(folderForImportantOutputFiles.getPath() + "/alignmentReport.txt");
        File fileMappingsHistogramPDF = new File(folderForImportantOutputFiles.getPath() + "/alignmentsByScore.histograms.pdf");
        File fileMappingsHistogramCumulativePDF = new File(folderForImportantOutputFiles.getPath() + "/alignmentsByScore.histograms.cumulative.pdf");

        // Split the read sequences
        File filesSplitReadsCSFasta[];
        int lengthsOfSplitReads[];
        int positionsStartOfSplitReadsInFullRead[];
        int countsMaxErrorsAllowedWhenMatching[];
        String masksByReadSplit[];
        int numberOfReads = -1;
        FastaDatabase fastaDBOfAllReads = new FastaDatabase(FILE_FULL_LENGTH_READS_CSFASTA);
        if (LENGTH_OF_READS >= 25) {   // split the reads into pieces
            if (LENGTH_OF_FIRST_PART_OF_READ > 0 && LENGTH_OF_LAST_PART_OF_READ > 0) {
                filesSplitReadsCSFasta = new File[2];
                filesSplitReadsCSFasta[0] = new File(folderForTempFiles.getPath() + "/reads.first_" + LENGTH_OF_FIRST_PART_OF_READ + EXTENSION_FOR_READS_FILE);
                filesSplitReadsCSFasta[1] = new File(folderForTempFiles.getPath() + "/reads.last_" + LENGTH_OF_LAST_PART_OF_READ + EXTENSION_FOR_READS_FILE);
                if (RUN_READ_SPLITTING) {
                    System.out.println("\nStarted read splitting at " + new Date(System.currentTimeMillis()));
                    System.out.println("\nSplitting the read sequences:\n" +
                            " -->First part of sequence will be colors 1 to " + LENGTH_OF_FIRST_PART_OF_READ + ":\t" + filesSplitReadsCSFasta[0].getPath() + "\n" +
                            " -->Last part of sequence will be colors " + (LENGTH_OF_READS - LENGTH_OF_LAST_PART_OF_READ +1) + " to " + LENGTH_OF_READS + ":\t" + filesSplitReadsCSFasta[1].getPath());
                    fastaDBOfAllReads.generateSplitVersionsToFiles(LENGTH_OF_FIRST_PART_OF_READ +1, LENGTH_OF_LAST_PART_OF_READ, "T0", filesSplitReadsCSFasta[0], filesSplitReadsCSFasta[1]);
                    numberOfReads = fastaDBOfAllReads.getNumberOfSequencesInDatabase();
                    writeNumberToFile(numberOfReads, fileReadsCount);
                    System.out.println("\nFinished read splitting at " + new Date(System.currentTimeMillis()));
                } else if (filesSplitReadsCSFasta[0].exists() && filesSplitReadsCSFasta[1].exists()) {
                    BufferedReader readerReadsCountFile = new BufferedReader(new FileReader(fileReadsCount));
                    numberOfReads = Integer.parseInt(readerReadsCountFile.readLine());
                    readerReadsCountFile.close();
                } else
                    throw new Exception("The split reads files do not exist.\nPlease re-run " + NAME_OF_TOOL + " with RUN_READ_SPLITTING turned on.");

                lengthsOfSplitReads = new int[]{LENGTH_OF_FIRST_PART_OF_READ, LENGTH_OF_LAST_PART_OF_READ +1};
                positionsStartOfSplitReadsInFullRead = new int[]{1, LENGTH_OF_READS - LENGTH_OF_LAST_PART_OF_READ}; // the position of the last part of a read is shifted down by one from the
                                                                                                                    // true position because we must add the faux "0" color at the front
                countsMaxErrorsAllowedWhenMatching = new int[]{MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING, MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING};

                masksByReadSplit = new String[]{ MASK_OF_READS.substring(0, LENGTH_OF_FIRST_PART_OF_READ),
                                                 "0" + MASK_OF_READS.substring(LENGTH_OF_READS - LENGTH_OF_LAST_PART_OF_READ)};
                System.out.println("\nmasksByReadSplit_0:\t" + masksByReadSplit[0]);
                System.out.println("masksByReadSplit_1:\t" + masksByReadSplit[1]);


            } else if (LENGTH_OF_LAST_PART_OF_READ <= 0) {

                filesSplitReadsCSFasta = new File[1];
                filesSplitReadsCSFasta[0] = new File(folderForTempFiles.getPath() + "/reads.first_" + LENGTH_OF_FIRST_PART_OF_READ + EXTENSION_FOR_READS_FILE);
                if (RUN_READ_SPLITTING) {
                    System.out.println("\nStarted read splitting at " + new Date(System.currentTimeMillis()));
                    System.out.println("\nSplitting the read sequences:\n" +
                            " -->First part of sequence will be colors 1 to " + LENGTH_OF_FIRST_PART_OF_READ + ":\t" + filesSplitReadsCSFasta[0].getPath() + "\n" +
                            " -->Last part of sequence will not be run.");
                    fastaDBOfAllReads.generateSplitVersionsToFiles(LENGTH_OF_FIRST_PART_OF_READ +1, LENGTH_OF_LAST_PART_OF_READ, "T0", filesSplitReadsCSFasta[0], null);
                    numberOfReads = fastaDBOfAllReads.getNumberOfSequencesInDatabase();
                    writeNumberToFile(numberOfReads, fileReadsCount);
                    System.out.println("\nFinished read splitting at " + new Date(System.currentTimeMillis()));
                } else if (filesSplitReadsCSFasta[0].exists()) {
                    BufferedReader readerReadsCountFile = new BufferedReader(new FileReader(fileReadsCount));
                    numberOfReads = Integer.parseInt(readerReadsCountFile.readLine());
                    readerReadsCountFile.close();
                } else
                    throw new Exception("The split reads files do not exist.\nPlease re-run " + NAME_OF_TOOL + " with " + RUN_READ_SPLITTING + " turned on.");

                lengthsOfSplitReads = new int[]{LENGTH_OF_FIRST_PART_OF_READ};
                positionsStartOfSplitReadsInFullRead = new int[]{1}; // the position of the last part of a read is shifted down by one from the
                                                                     // true position because we must add the faux "0" color at the front
                countsMaxErrorsAllowedWhenMatching = new int[]{MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING};

                masksByReadSplit = new String[]{ MASK_OF_READS.substring(0, LENGTH_OF_FIRST_PART_OF_READ)};
                System.out.println("\nmasksByReadSplit_0:\t" + masksByReadSplit[0]);

            } else
                throw new Exception("First parts of reads are too short (L = " + LENGTH_OF_FIRST_PART_OF_READ + ").");

        } else
            throw new Exception("Reads are too short (L = " + LENGTH_OF_READS + ").");

        JobSubmissionParameters parameterTemplate = new JobSubmissionParameters();
        parameterTemplate.setEnvironment(SCHEDULING_ENVIRONMENT);
        parameterTemplate.setQueueName(NAME_OF_QUEUE);
        parameterTemplate.setRerunnable(JOBS_RERUNNABLE);
        parameterTemplate.setResourceString(SCHEDULER_RESOURCE_REQUIREMENTS);
        parameterTemplate.setAdditionalOptions(ADDITIONAL_SCHEDULER_OPTIONS);

        
        
        // Filter reads that match to a database of "junk" sequences (e.g., rRNA, primers, repeats) prior to mapping to the reference genome
        if (RUN_READ_FILTERING) {

            System.out.println("\nStarted read filtering at " + new Date(System.currentTimeMillis()));
            System.out.println("\nReads will be filtered against the following fasta file:\t" + FILE_FILTER_SEQUENCES_FASTA.getPath());

            ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();
            
            long memoryRequiredForMapReads = (long)Math.ceil(ReadMapper.calculateMemoryRequiredByMapReadsInBytes(new FastaDatabase(FILE_FILTER_SEQUENCES_FASTA).getTotalLengthOfSequence()) * MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR);
            
            JobSubmissionParameters params = parameterTemplate.clone();
            params.setMemoryRequirement(memoryRequiredForMapReads);
            ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);

            File[] filesMAByReadSequenceSplit = new File[filesSplitReadsCSFasta.length];
            File[] filesFilteredReadsByReadSequenceSplit = new File[filesSplitReadsCSFasta.length];

            for (int indexSplitReadsFile = 0; indexSplitReadsFile < filesSplitReadsCSFasta.length; indexSplitReadsFile++) {

                File folderForFilterIO = new File(folderForFilterMapping.getPath() + "/reads_" + indexSplitReadsFile + "/ref_0/");
                if (!folderForFilterIO.exists())
                    folderForFilterIO.mkdirs();

                ReadMapper mapper = new ReadMapper(fileMapReadsExe, FILE_FILTER_SEQUENCES_FASTA, filesSplitReadsCSFasta[indexSplitReadsFile],
                                                   lengthsOfSplitReads[indexSplitReadsFile], masksByReadSplit[indexSplitReadsFile]);


                int numberOfReadsPerJob = (int)Math.ceil((double)numberOfReads / NUMBER_OF_READS_FILE_PARTITIONS);
                for (int indexOfRead = 0; indexOfRead < numberOfReads; indexOfRead += numberOfReadsPerJob) {

                    int maxReadIndexInCurrentRange = Math.min(indexOfRead + numberOfReadsPerJob, numberOfReads) -1;

                    String nameOfMAFile = filesSplitReadsCSFasta[indexSplitReadsFile].getName()
                                            + "." + lengthsOfSplitReads[indexSplitReadsFile] + "." + countsMaxErrorsAllowedWhenMatching[indexSplitReadsFile];
                    if (VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION)
                        nameOfMAFile += ".adj_valid";
                    nameOfMAFile += "." + indexOfRead + "_" + maxReadIndexInCurrentRange + EXTENSION_FOR_MAPREADS_OUTPUT_FILE;

                    filesMAByReadSequenceSplit[indexSplitReadsFile] = new File(folderForFilterIO.getPath() + "/" + nameOfMAFile);

                    mapper.startMapReadsJob(true, indexOfRead, maxReadIndexInCurrentRange,
                                            MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING, false, true, 1,
                                            clusterInterface, "wt_filter.",
                                            FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES, FOLDER_CONTAINING_SCHEMA_FILES, folderForFilterIO, filesMAByReadSequenceSplit[indexSplitReadsFile], true);
                     
                    listOfFilesThatAreNoLongerNeeded.add(filesMAByReadSequenceSplit[indexSplitReadsFile]);
                }

                filesFilteredReadsByReadSequenceSplit[indexSplitReadsFile] = new File(filesSplitReadsCSFasta[indexSplitReadsFile].getPath().replaceAll(EXTENSION_FOR_READS_FILE + "$", EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE));
            }

            clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterFilteringJobSubmissionScript);
            clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterFilteringJobRemovalScript);
            clusterInterface.writeMasterListOfJobOutputFilesFromLog(fileMasterFilteringListOfScriptOutputFiles);

            fileMasterFilteringJobSubmissionScript.setExecutable(true, true);
            fileMasterFilteringJobRemovalScript.setExecutable(true, true);

            // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
            System.out.print("\nWaiting for filtering jobs to finish...");
            while (!clusterInterface.checkIfLoggedJobsComplete())
                Thread.sleep(5000);

            System.out.println("OK");

            ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
            if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
                System.out.println("The following jobs failed:");
                Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
                while (iteratorOverFailFiles.hasNext())
                    System.out.println(((File) iteratorOverFailFiles.next()).getPath());

                throw new Exception("One or more filtering jobs failed.");
            }

            System.out.println("\nTagging filtered reads and generating a filter report.");

            FilteringTask.processFilteredReadsAndWriteReport(true,
                                                MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING,
                                                null,
                                                FILE_FILTER_SEQUENCES_FASTA,
                                                filesMAByReadSequenceSplit,
                                                filesFilteredReadsByReadSequenceSplit,
                                                fileFilterReportTxt);

            if (FILTERING_MODE != FilteringMode.OFF)
                for (int indexSplitReadsFile = 0; indexSplitReadsFile < filesSplitReadsCSFasta.length; indexSplitReadsFile++)
                    filesSplitReadsCSFasta[indexSplitReadsFile] = filesFilteredReadsByReadSequenceSplit[indexSplitReadsFile];


            if (DELETE_INTERMEDIATE_FILES || COMPRESS_INTERMEDIATE_FILES) {
                if (COMPRESS_INTERMEDIATE_FILES && DELETE_INTERMEDIATE_FILES)  System.out.print("\nCompressing and deleting intermediate files...");
                else if (COMPRESS_INTERMEDIATE_FILES)   System.out.print("\nCompressing intermediate files...");
                else    System.out.print("\nDeleting intermediate files...");
                Iterator<File> iteratorFilesNoLongerNeeded = listOfFilesThatAreNoLongerNeeded.iterator();
                while (iteratorFilesNoLongerNeeded.hasNext()) {
                    File file = iteratorFilesNoLongerNeeded.next();
                    if (COMPRESS_INTERMEDIATE_FILES)
                        CompressionUtilities.compressFile(file, new File(file.getPath() + EXTENSION_FOR_COMPRESSED_FILE));
                    if (DELETE_INTERMEDIATE_FILES)
                        file.delete();
                }
                System.out.println("OK");
            }

            System.out.println("\nFinished read filtering at " + new Date(System.currentTimeMillis()));

        } else if (FILTERING_MODE != FilteringMode.OFF) {    // if using pre-existing filtered reads diles, check that all the filter mapping jobs completed successfully and that the merge exists

            System.out.print("\nChecking that filtering jobs from previous run completed successfully...");
            checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(fileMasterFilteringListOfScriptOutputFiles, "\t", 0));
            System.out.println(" OK");

            for (int indexReadsFile = 0; indexReadsFile < filesSplitReadsCSFasta.length; indexReadsFile++) {
                System.out.println("Checking that post-filter reads file exist: " + filesSplitReadsCSFasta[indexReadsFile]);

                File fileMergedFilteredOutput = new File(filesSplitReadsCSFasta[indexReadsFile].getPath().replaceAll(EXTENSION_FOR_READS_FILE + "$", EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE));
                if (!fileMergedFilteredOutput.exists())
                    throw new Exception("Could not find filtered reads file " + fileMergedFilteredOutput + ".  Please re-run " + NAME_OF_TOOL + " with RUN_READ_FILTERING on.");

                filesSplitReadsCSFasta[indexReadsFile] = fileMergedFilteredOutput;
            }

        } else {

            System.out.println("FILTER_READS_WITH_ALL_SPLITS_MAPPING_TO_THE_FILTER_SEQUENCES and FILTER_READS_WITH_ONE_OR_MORE_SPLITS_MAPPING_TO_THE_FILTER_SEQUENCES are off, therefore will not filter reads.");

        }



        // Split the reference file into pieces
        File filesReferenceFastaPartitions[];
        List<Integer> listsOfFullRefSequenceNumbersOrderedForEachReferencePartition[];   // tracks the reference header ordering, for aggregating data later
        if (RUN_REFERENCE_PARTITIONING) {

            System.out.println("\nStarted reference partitioning at " + new Date(System.currentTimeMillis()));
            System.out.println("\nReference file will be partitioned into files with no more than " + (new DecimalFormat("0.00E0")).format(MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES) + " bases each.");
            
            FastaDatabase fastaDBOfFullReference = new FastaDatabase(FILE_REFERENCE_FASTA);
            @SuppressWarnings("unchecked")
            ArrayList<HashSet<String>> listOfReferenceHeaderSets = fastaDBOfFullReference.getSetsOfSequenceHeadersEachTotalingLessThanNMegabases(MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES);

            filesReferenceFastaPartitions = new File[listOfReferenceHeaderSets.size()];
            for (int i = 0; i < filesReferenceFastaPartitions.length; i++)
                filesReferenceFastaPartitions[i] = new File(folderForTempFiles.getPath() + "/reference." + i + ".fasta");

            listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = FastaDatabase.partitionFastaFileToFastaFilesByHeaderSets(fileReferencePartitionTable, filesReferenceFastaPartitions, fastaDBOfFullReference, listOfReferenceHeaderSets);

            System.out.println("\nFinished reference partitioning at " + new Date(System.currentTimeMillis()));

        } else {

            if (!fileReferencePartitionTable.exists())
                throw new Exception("The reference partition table file must exist: " + fileReferencePartitionTable.getPath() + "\nPlease re-run with RUN_REFERENCE_PARTITIONING turned on.");

            listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = loadPartitionTable(fileReferencePartitionTable);

            int maxRefPartitionIndex = getMaxPartitionIndexFromPartitionTableFile(fileReferencePartitionTable);
            filesReferenceFastaPartitions = new File[maxRefPartitionIndex +1];
            for (int i = 0; i < filesReferenceFastaPartitions.length; i++)
                filesReferenceFastaPartitions[i] = new File(folderForTempFiles.getPath() + "/reference." + i + ".fasta");

        }



        if (RUN_MAPPING) {

            System.out.println("\nStarted mapping jobs at " + new Date(System.currentTimeMillis()));

            
            long memoryRequiredForMapReads = (long)Math.ceil(ReadMapper.calculateMemoryRequiredByMapReadsInBytes(MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES)*MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR);
            JobSubmissionParameters params = parameterTemplate.clone();
            params.setMemoryRequirement(memoryRequiredForMapReads);
            ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);

            for (int indexReadSplit = 0; indexReadSplit < filesSplitReadsCSFasta.length; indexReadSplit++) {

                for (int indexRefFile = 0; indexRefFile < filesReferenceFastaPartitions.length; indexRefFile++) {

                    File subfolderForOutputFiles = new File(folderForMappingAndExtension.getPath() + "/reads_" + indexReadSplit + "/ref_" + indexRefFile);
                    if (!subfolderForOutputFiles.exists())
                        subfolderForOutputFiles.mkdirs();

                    ReadMapper mapper = new ReadMapper(fileMapReadsExe, filesReferenceFastaPartitions[indexRefFile], filesSplitReadsCSFasta[indexReadSplit],
                                                       lengthsOfSplitReads[indexReadSplit], masksByReadSplit[indexReadSplit]);

                    int numberOfReadsPerJob = (int)Math.ceil((double)numberOfReads / NUMBER_OF_READS_FILE_PARTITIONS);
                    for (int indexOfRead = 0; indexOfRead < numberOfReads; indexOfRead += numberOfReadsPerJob) {

                        int maxReadIndexInCurrentRange = Math.min(indexOfRead + numberOfReadsPerJob, numberOfReads) -1;

                        String nameOfMAFile = filesSplitReadsCSFasta[indexReadSplit].getName()
                                                + "." + lengthsOfSplitReads[indexReadSplit] + "." + countsMaxErrorsAllowedWhenMatching[indexReadSplit];
                        if (VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION)
                            nameOfMAFile += ".adj_valid";
                        nameOfMAFile += "." + indexOfRead + "_" + maxReadIndexInCurrentRange + EXTENSION_FOR_MAPREADS_OUTPUT_FILE;

                        File fileMA = new File(subfolderForOutputFiles.getPath() + "/" + nameOfMAFile);

                        mapper.startMapReadsJob(true, indexOfRead, maxReadIndexInCurrentRange,
                                                countsMaxErrorsAllowedWhenMatching[indexReadSplit], VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION, MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION,
                                                (MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ +1),
                                                clusterInterface, "wt_map.",
                                                FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES, FOLDER_CONTAINING_SCHEMA_FILES, subfolderForOutputFiles, fileMA, true);
                    }


                }

            }

            clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterMappingJobSubmissionScript);
            clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterMappingJobRemovalScript);
            clusterInterface.writeMasterListOfJobOutputFilesFromLog(fileMasterMappingListOfScriptOutputFiles);

            fileMasterMappingJobSubmissionScript.setExecutable(true, true);
            fileMasterMappingJobRemovalScript.setExecutable(true, true);

            // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
            System.out.print("\nWaiting for mapping jobs to finish...");
            while (!clusterInterface.checkIfLoggedJobsComplete())
                Thread.sleep(5000);

            System.out.println("OK");

            ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
            if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
                System.out.println("The following jobs failed:");
                Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
                while (iteratorOverFailFiles.hasNext())
                    System.out.println(( iteratorOverFailFiles.next()).getPath());

                throw new Exception("One or more mapping jobs failed.");
            }

            System.out.println("\nFinished mapping jobs at " + new Date(System.currentTimeMillis()));

        } else {    // if skipping mapping, check that all the mapping jobs completed successfully

            System.out.print("\nChecking that mapping jobs from previous run completed successfully...");
            checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(fileMasterMappingListOfScriptOutputFiles, "\t", 0));
            System.out.println(" OK");
        }




        if (RUN_EXTENSION) {

            System.out.println("\nStarted extension jobs at " + new Date(System.currentTimeMillis()));

            ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();

            long memoryRequiredByExtension = (long)Math.ceil(ReadMappingExtender.calculateMemoryRequiredByExtendMapReadsInBytes(MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES) * MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR);

            JobSubmissionParameters params = parameterTemplate.clone();
            params.setMemoryRequirement(memoryRequiredByExtension);
            ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);

            for (int indexReadsFile = 0; indexReadsFile < filesSplitReadsCSFasta.length; indexReadsFile++) {

                for (int indexRefFile = 0; indexRefFile < filesReferenceFastaPartitions.length; indexRefFile++) {

                    File subfolderForOutputFiles = new File(folderForMappingAndExtension.getPath() + "/reads_" + indexReadsFile + "/ref_" + indexRefFile);
                    if (!subfolderForOutputFiles.exists())
                        throw new Exception("Could not find folder with MA output files: " + subfolderForOutputFiles.getPath() + "\nPlease ensure that mapping was run prior to extension.");

                    int numberOfReadsPerJob = (int)Math.ceil((double)numberOfReads / NUMBER_OF_READS_FILE_PARTITIONS);
                    for (int indexOfRead = 0; indexOfRead < numberOfReads; indexOfRead += numberOfReadsPerJob) {

                        int maxReadIndexInCurrentRange = Math.min(indexOfRead + numberOfReadsPerJob, numberOfReads) -1;

                        String prefixOfMAFile = filesSplitReadsCSFasta[indexReadsFile].getName()
                                                + "." + lengthsOfSplitReads[indexReadsFile] + "." + countsMaxErrorsAllowedWhenMatching[indexReadsFile];
                        if (VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION)
                            prefixOfMAFile += ".adj_valid";
                        String nameOfMAFile  = prefixOfMAFile + "." + indexOfRead + "_" + maxReadIndexInCurrentRange + EXTENSION_FOR_MAPREADS_OUTPUT_FILE;

                        File fileMA = new File(subfolderForOutputFiles.getPath() + "/" + nameOfMAFile);
                        File fileMAX = new File(subfolderForOutputFiles.getPath() + "/" + nameOfMAFile.replaceFirst(EXTENSION_FOR_MAPREADS_OUTPUT_FILE + "$", EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE));

                        if (!fileMA.exists()) {
                            File fileMACompressed = new File(fileMA.getPath() + EXTENSION_FOR_COMPRESSED_FILE);
                            if (fileMACompressed.exists()) {
                                System.out.print("\nDecompressing intermediate files...");
                                CompressionUtilities.decompressFiles(fileMACompressed, fileMACompressed.getParentFile());
                                System.out.println("OK");
                            } else
                                throw new Exception("Could not find MA file:\t" + fileMA.getPath() + "\nPlease re-run " + NAME_OF_TOOL + " with RUN_MAPPING turned on.");
                        }

                        // remove merged MAX files if found because these should be re-created from the new MAX files
                        File fileMergedMax = new File(subfolderForOutputFiles.getPath() + "/" + prefixOfMAFile + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
                        if (fileMergedMax.exists()) fileMergedMax.delete();


                        ReadMappingExtender extender = new ReadMappingExtender(fileExtendMapReadsExe, filesReferenceFastaPartitions[indexRefFile],
                                                                               fileMA, FILE_FULL_LENGTH_READS_CSFASTA, LENGTH_OF_READS);

                        extender.startExtendMappedReadsJob(SCORE_OF_MATCH_FOR_EXTENSION, SCORE_OF_MISMATCH_FOR_EXTENSION,
                                                            VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION,
                                                            MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION,
                                                            positionsStartOfSplitReadsInFullRead[indexReadsFile],
                                                            clusterInterface, "wt_extension.",
                                                            FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES, subfolderForOutputFiles,
                                                            fileMAX, new File(fileMAX.getPath() + EXTENSION_FOR_INDEX_FILE),
                                                            true, MASK_OF_READS);

                        listOfFilesThatAreNoLongerNeeded.add(fileMA);
                    }

                }

            }

            clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterExtensionJobSubmissionScript);
            clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterExtensionJobRemovalScript);
            clusterInterface.writeMasterListOfJobOutputFilesFromLog(fileMasterExtensionListOfScriptOutputFiles);

            fileMasterExtensionJobSubmissionScript.setExecutable(true, true);
            fileMasterExtensionJobRemovalScript.setExecutable(true, true);

            // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
            System.out.print("\nWaiting for extension jobs to finish...");
            while (!clusterInterface.checkIfLoggedJobsComplete())
                Thread.sleep(5000);

            System.out.println("OK");

            ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
            if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
                System.out.println("The following jobs failed:");
                Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
                while (iteratorOverFailFiles.hasNext())
                    System.out.println(( iteratorOverFailFiles.next()).getPath());

                throw new Exception("One or more extension jobs failed.");
            }


            if (DELETE_INTERMEDIATE_FILES || COMPRESS_INTERMEDIATE_FILES) {
                if (COMPRESS_INTERMEDIATE_FILES && DELETE_INTERMEDIATE_FILES)  System.out.print("\nCompressing and deleting intermediate files...");
                else if (COMPRESS_INTERMEDIATE_FILES)   System.out.print("\nCompressing intermediate files...");
                else    System.out.print("\nDeleting intermediate files...");
                Iterator<File> iteratorFilesNoLongerNeeded = listOfFilesThatAreNoLongerNeeded.iterator();
                while (iteratorFilesNoLongerNeeded.hasNext()) {
                    File file = iteratorFilesNoLongerNeeded.next();
                    if (COMPRESS_INTERMEDIATE_FILES)
                        CompressionUtilities.compressFile(file, new File(file.getPath() + EXTENSION_FOR_COMPRESSED_FILE));
                    if (DELETE_INTERMEDIATE_FILES)
                        file.delete();
                }
                System.out.println("OK");
            }


            System.out.println("\nFinished extension jobs at " + new Date(System.currentTimeMillis()));

        } else {    // if skipping extension, check that all the mapping jobs completed successfully

            System.out.print("\nChecking that extension jobs from previous run completed successfully...");
            checkThatAllJobsCompletedSuccessfully(TextFileUtilities.loadSetFromFile(fileMasterExtensionListOfScriptOutputFiles, "\t", 0));
            System.out.println(" OK");
        }



        if (RUN_MERGE) {

            System.out.println("\nStarted merge of MAX files at " + new Date(System.currentTimeMillis()));

            ArrayList<File> listOfFilesThatAreNoLongerNeeded = new ArrayList<File>();

            int maxMinPossibleAlignmentScore = Integer.MIN_VALUE;
            for (int indexReadSplit = 0; indexReadSplit < lengthsOfSplitReads.length; indexReadSplit++) {
                int numberOfUnmaskedPositions = masksByReadSplit[indexReadSplit].replaceAll("0", "").length();
                maxMinPossibleAlignmentScore = Math.max(maxMinPossibleAlignmentScore,
                                                (numberOfUnmaskedPositions - countsMaxErrorsAllowedWhenMatching[indexReadSplit]) * SCORE_OF_MATCH_FOR_EXTENSION + countsMaxErrorsAllowedWhenMatching[indexReadSplit] * SCORE_OF_MISMATCH_FOR_EXTENSION );
            }
            System.out.println("\nmaxMinPossibleAlignmentScore:\t" + maxMinPossibleAlignmentScore);

            // locate all the MAX files and merge them at the reads file partition level, leaving reference file partitions and read sequence splits intact
            File[][] filesMAXByRefPartitionByReadSequenceSplit = new File[filesReferenceFastaPartitions.length][filesSplitReadsCSFasta.length];
            File[][] filesMAXIndexByRefPartitionByReadSequenceSplit = new File[filesReferenceFastaPartitions.length][filesSplitReadsCSFasta.length];
            for (int indexReferenceFilePartition = 0; indexReferenceFilePartition < filesMAXByRefPartitionByReadSequenceSplit.length; indexReferenceFilePartition++)
                for (int indexReadsFilePartition = 0; indexReadsFilePartition < filesMAXByRefPartitionByReadSequenceSplit[0].length; indexReadsFilePartition++) {

                    File folderCurrentPartition = new File(folderForMappingAndExtension.getPath() + "/reads_" + indexReadsFilePartition + "/ref_" + indexReferenceFilePartition + "/");

                    // does the MAX file merged across partitions already exist?
                    File filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[] = folderCurrentPartition.listFiles(FILENAME_FILTER_MERGED_MAX);
                    if (filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged.length == 1)
                        filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[0];
                    else {
                        // don't know why there would be more than one file, but if so cleanup
                        for (int i = 0; i < filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged.length; i++)
                            filesMAXForThisRefPartitionAndThisReadSequenceSplitMerged[i].delete();

                        // if not, has it been compressed?
                        File filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[] = folderCurrentPartition.listFiles(FILENAME_FILTER_MERGED_MAX_COMPRESSED);
                        if (filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped.length == 1) {
                            System.out.print("\nDecompressing intermediate files...");
                            CompressionUtilities.decompressFiles(filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[0], folderCurrentPartition);
                            filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = folderCurrentPartition.listFiles(FILENAME_FILTER_MERGED_MAX)[0];
                            System.out.println("OK");
                        } else {

                            // don't know why there would be more than one file, but if so cleanup
                            for (int i = 0; i < filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped.length; i++)
                                filesMAXForThisRefPartitionAndThisReadSequenceSplitMergedZipped[i].delete();

                            // if not compressed, are the unmerged MAX files available to make the merge anew?
                            File filesMAXForThisRefPartitionAndThisReadSequenceSplit[] = folderCurrentPartition.listFiles(FILENAME_FILTER_MAX);
                            if (filesMAXForThisRefPartitionAndThisReadSequenceSplit.length == 1) {
                                int indexOfSecondToLastPeriod = filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().lastIndexOf(".", filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().length() - EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE.length() -1);
                                filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath().substring(0, indexOfSecondToLastPeriod) + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE);
                                filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition] = new File(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition].getPath() + EXTENSION_FOR_INDEX_FILE);
                                filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].renameTo(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);
                                new File(filesMAXForThisRefPartitionAndThisReadSequenceSplit[0].getPath() + EXTENSION_FOR_INDEX_FILE).renameTo(filesMAXIndexByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);

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
                                throw new Exception("Could not find MAX files in:\t" + folderCurrentPartition.getPath() + "\nPlease re-run " + NAME_OF_TOOL + " with RUN_EXTENSION turned on.");

                        }
                    }
                    listOfFilesThatAreNoLongerNeeded.add(filesMAXByRefPartitionByReadSequenceSplit[indexReferenceFilePartition][indexReadsFilePartition]);
                }

            JobSubmissionParameters params = parameterTemplate.clone();
            params.setMemoryRequirement(MAX_MEMORY_PER_JOB_IN_BYTES.longValue());
            ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
            // Need to break up the reads into a new set of partitions that can fit into memory on each node (for sorting).  This is a worst case guess for now.
            // FIXME: IS the calculation robust?
            int numberOfReadsPerPartition = (int) (MAX_MEMORY_PER_JOB_IN_BYTES/MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR / 4000.0 *  10.0 / MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ);   //1500.0);
            if (NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE != null)
                numberOfReadsPerPartition = NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE;
            File filesOutputForEachReadsFilePartition[] = new File[(int)Math.ceil((double)numberOfReads / numberOfReadsPerPartition)];
            File filesSortedOutputForEachReadsFilePartition[] = new File[filesOutputForEachReadsFilePartition.length];
            File filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[] = new File[filesOutputForEachReadsFilePartition.length];
            File filesSerializedReportOfExtendedReadMappingsOut[] = new File[filesOutputForEachReadsFilePartition.length];
            for (int indexOfPartition = 0; indexOfPartition < filesOutputForEachReadsFilePartition.length; indexOfPartition++) {

                int indexOfFirstReadToMerge = indexOfPartition * numberOfReadsPerPartition;
                int indexOfLastReadToMerge = Math.min(indexOfFirstReadToMerge + numberOfReadsPerPartition -1, numberOfReads -1);

                System.out.println("\n -->Merging MAX files by reference partition and read sequence split for reads " + indexOfFirstReadToMerge + " to " + indexOfLastReadToMerge);

                File folderForMergingThisReadPartition = new File(folderForMerge.getPath() + "/reads_" + indexOfFirstReadToMerge + "_to_" + indexOfLastReadToMerge);
                if (!folderForMergingThisReadPartition.exists())
                    folderForMergingThisReadPartition.mkdirs();

                filesOutputForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE.substring(1) + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE);
                listOfFilesThatAreNoLongerNeeded.add(filesOutputForEachReadsFilePartition[indexOfPartition]);

                filesSortedOutputForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + EXTENSION_FOR_SORTED_OUTPUT_FILE.substring(1) + EXTENSION_FOR_MAPREADS_PLUS_EXTENSION_OUTPUT_FILE + EXTENSION_FOR_MERGED_OUTPUT_FILE + EXTENSION_FOR_FILTERED_OUTPUT_FILE + EXTENSION_FOR_READS_FILE);
                listOfFilesThatAreNoLongerNeeded.add(filesSortedOutputForEachReadsFilePartition[indexOfPartition]);

                filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/" + EXTENSION_FOR_READS_WITHOUT_MIN_SCORING_ALIGNMENT_FILE.substring(1) + EXTENSION_FOR_READS_FILE);
                listOfFilesThatAreNoLongerNeeded.add(filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition]);

                filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition] = new File(folderForMergingThisReadPartition.getPath() + "/reportOfExtendedReadMappingsOut." + System.currentTimeMillis() + ".obj");
                listOfFilesThatAreNoLongerNeeded.add(filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition]);

                ExtendedReadMappingsMerger.startMergeExtendedReadMappingsJob(MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT,
                                                                            MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS,
                                                                            MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ,
                                                                            MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ,
                                                                            maxMinPossibleAlignmentScore,
                                                                            LENGTH_OF_READS * SCORE_OF_MATCH_FOR_EXTENSION,
                                                                            indexOfFirstReadToMerge,
                                                                            indexOfLastReadToMerge,
                                                                            0,1,
                                                                            FILTERING_MODE, 
                                                                            EXTENSION_FOR_INDEX_FILE,
                                                                            listsOfFullRefSequenceNumbersOrderedForEachReferencePartition,
                                                                            filesMAXByRefPartitionByReadSequenceSplit,
                                                                            null,
                                                                            //reportOfExtendedReadMappings,
                                                                            clusterInterface,
                                                                            "wt_merge.",
                                                                            folderContainingJavaLibraries.getPath() + "/'*':"+WHOLE_TRANSCRIPTOME_JAR.getPath(),
                                                                            FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES,
                                                                            folderForMergingThisReadPartition,
                                                                            folderForMergingThisReadPartition,
                                                                            filesSerializedReportOfExtendedReadMappingsOut[indexOfPartition],
                                                                            filesOutputForEachReadsFilePartition[indexOfPartition],
                                                                            filesSortedOutputForEachReadsFilePartition[indexOfPartition],
                                                                            filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition[indexOfPartition],
                                                                            true);


            }


            clusterInterface.writeMasterJobSubmissionFileFromLog(fileMasterMergeJobSubmissionScript);
            clusterInterface.writeMasterJobRemovalFileFromLog(fileMasterMergeJobRemovalScript);
            clusterInterface.writeMasterListOfJobOutputFilesFromLog(fileMasterMergeListOfScriptOutputFiles);

            fileMasterMergeJobSubmissionScript.setExecutable(true, true);
            fileMasterMergeJobRemovalScript.setExecutable(true, true);

            // FIXME: If one job fails before all jobs complete it would be beneficial to report this earlier.
            System.out.print("\nWaiting for merge jobs to finish...");
            while (!clusterInterface.checkIfLoggedJobsComplete())
                Thread.sleep(5000);

            System.out.println("OK");

            ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
            if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
                System.out.println("The following jobs failed:");
                Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
                while (iteratorOverFailFiles.hasNext())
                    System.out.println(( iteratorOverFailFiles.next()).getPath());

                throw new Exception("One or more merge jobs failed.");
            }
            

            System.out.println("\n -->Merging unsorted MAX files, by reads file partition.");
            Utilities.concatenateFiles(filesOutputForEachReadsFilePartition, fileFullyMergedMAXcsfasta, false);
            ExtendedReadMappingsReport reportOfExtendedReadMappings = ExtendedReadMappingsReport.combineSerializedMappingReports(filesSerializedReportOfExtendedReadMappingsOut);


            System.out.println("\n -->Merging sorted MAX files, by reads file partition.");
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileFullyMergedAndSortedMAXcsfasta));
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


            System.out.println("\n -->Merging reads w/o min scoring alignments files, by reads file partition.");
            Utilities.concatenateFiles(filesReadsWithoutMinScoringAlignmentForEachReadsFilePartition, fileReadsWithoutMinScoringAlignmentCSfasta, false);


            if (DELETE_INTERMEDIATE_FILES || COMPRESS_INTERMEDIATE_FILES) {
                if (COMPRESS_INTERMEDIATE_FILES && DELETE_INTERMEDIATE_FILES)  System.out.print("\nCompressing and deleting intermediate files...");
                else if (COMPRESS_INTERMEDIATE_FILES)   System.out.print("\nCompressing intermediate files...");
                else    System.out.print("\nDeleting intermediate files...");
                Iterator<File> iteratorFilesNoLongerNeeded = listOfFilesThatAreNoLongerNeeded.iterator();
                while (iteratorFilesNoLongerNeeded.hasNext()) {
                    File file = iteratorFilesNoLongerNeeded.next();
                    if (COMPRESS_INTERMEDIATE_FILES)
                        CompressionUtilities.compressFile(file, new File(file.getPath() + EXTENSION_FOR_COMPRESSED_FILE));
                    if (DELETE_INTERMEDIATE_FILES)
                        file.delete();
                }
                System.out.println("OK");
            }


            System.out.println("\nFinished merge of MAX files at " + new Date(System.currentTimeMillis()) + "\n");

            try {
                reportOfExtendedReadMappings.getHistogramOfMappingsByScore().setReportRelativeFrequencies(true);
                reportOfExtendedReadMappings.getHistogramOfMappingsByScore().toPDF(fileMappingsHistogramPDF, false, false, false, null);
                reportOfExtendedReadMappings.getHistogramOfMappingsByScore().toPDF(fileMappingsHistogramCumulativePDF, false, true, false, null);

            } catch (Exception e) {
                System.out.println("Caught exception while trying to generate PDF.");
                e.printStackTrace();
            }

            reportOfExtendedReadMappings.writeReportToTextFile(fileMappingReportTxt);
            System.out.println("\n\n" + reportOfExtendedReadMappings.getReport1_1());

            
        }




    }

    private static void writeNumberToFile(int number, File file) throws IOException {
        BufferedWriter writerReadsCountFile = new BufferedWriter(new FileWriter(file));
        writerReadsCountFile.write(number + "");
        writerReadsCountFile.newLine();
        writerReadsCountFile.close();
    }

    private static void loadResourceBundle() {
    	resourceBundle = (PropertyResourceBundle) ResourceBundle.getBundle(WholeTranscriptomeAnalyzer.class.getPackage().getName()+".application");
    }

    private static void loadConfigurationFromFile(File fileConfiguration) throws Exception {

        System.out.println("Loading configuration file...\n");

        // SOME DEFAULTS THAT CAN NOT BE CHANGED
        //NUMBER_OF_NODES_ALLOCATED_PER_JOB = 1;
        JOBS_RERUNNABLE = false;

        // SOME DEFAULTS THAT CAN BE CHANGED
        //double maxMemoryPerJobInBytes = 4E9;
        //MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES = Math.min(1E9, maxMemoryPerJobInBytes / 6);
        NUMBER_OF_READS_FILE_PARTITIONS = 1;     // if more processors are available then split more times... just beware of memory requirements
                                                 // how can we automatically set this one?
        //NUMBER_OF_PROCESSORS_ALLOCATED_PER_NODE = 1;

        RUN_READ_SPLITTING = true;
        RUN_READ_FILTERING = false;
        RUN_REFERENCE_PARTITIONING = true;
        RUN_MAPPING = true;
        RUN_EXTENSION = true;
        RUN_MERGE = true;

        FILTERING_MODE = FilteringMode.OFF;

        COMPRESS_INTERMEDIATE_FILES = false;
        DELETE_INTERMEDIATE_FILES = true;

        MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ = 1;
        MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ = 10;
        VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION = true;
        MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION = false;

        SCORE_OF_MATCH_FOR_EXTENSION = 1;
        SCORE_OF_MISMATCH_FOR_EXTENSION = -1;

        MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING = 2;
        
        MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR = 1.0f;
        ADDITIONAL_SCHEDULER_OPTIONS = "";
        
        // Now read the configuration file.
        try {
            BufferedReader readerConfigurationFile = new BufferedReader(new FileReader(fileConfiguration));
            String line = null;
            while ((line = readerConfigurationFile.readLine()) != null) {

                if (line.length() > 0 && !line.startsWith("#")) {
                    String tokens[] = line.split("\t");
                    if (tokens.length >= 2) {

                        String nameOfParameter = tokens[0];
                        String valueOfParameter = tokens[1].trim();
                        if (nameOfParameter.equalsIgnoreCase("NAME_OF_QUEUE")) {
                            NAME_OF_QUEUE = valueOfParameter;
                        } else if (nameOfParameter.equalsIgnoreCase("SCHEDULING_ENVIRONMENT")) {
                        	SCHEDULING_ENVIRONMENT = valueOfParameter;
                        } else if (nameOfParameter.equalsIgnoreCase("NAME_OF_PROJECT")) {
                        	throw new Exception("The parameter: 'NAME_OF_PROJECT' is no longer valid.  Please remove it from your config file. Set project name using 'ADDITIONAL_SCHEDULER_OPTIONS'.");
                        } else if (nameOfParameter.equalsIgnoreCase("MAX_MEMORY_PER_JOB_IN_BYTES")) {

                            try {
                                MAX_MEMORY_PER_JOB_IN_BYTES = (long)Double.parseDouble(valueOfParameter);
                                MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES = ReadMapper.calculateMaxReferenceSizeForMapReads(MAX_MEMORY_PER_JOB_IN_BYTES);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MAX_MEMORY_PER_JOB_IN_BYTES must be numeric.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE")) {

                            try {
                                NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("NUMBER_OF_READS_FILE_PARTITIONS")) {

                        	throw new Exception("The parameter: 'NUMBER_OF_READS_FILE_PARTITIONS' is no longer valid.  Please remove it from your config file.");

//                            System.out.println("NUMBER_OF_READS_FILE_PARTITIONS can currently only be set to 1.");
//                            /*try {
//                                NUMBER_OF_READS_FILE_PARTITIONS = Integer.parseInt(valueOfParameter);
//                            } catch (NumberFormatException nfe) {
//                                System.out.println("NUMBER_OF_READS_FILE_PARTITIONS must be an integer.");
//                                throw nfe;
//                            } */

                        } else if (nameOfParameter.equalsIgnoreCase("NUMBER_OF_PROCESSORS_ALLOCATED_PER_NODE")) {

                        	throw new Exception("The parameter: 'NUMBER_OF_PROCESSORS_ALLOCATED_PER_NODE' is no longer valid.  Please remove it from your config file.  Please use 'SCHEDULER_RESOURCE_REQUIREMENTS' or 'ADDITIONAL_SCHEDULER_OPTIONS' to specify this setting for your scheduler.");

//                            try {
//                                NUMBER_OF_PROCESSORS_ALLOCATED_PER_NODE = Integer.parseInt(valueOfParameter);
//                            } catch (NumberFormatException nfe) {
//                                System.out.println("NUMBER_OF_PROCESSORS_ALLOCATED_PER_NODE must be an integer.");
//                                throw nfe;
//                            }

                        } else if (nameOfParameter.equalsIgnoreCase("RUN_READ_SPLITTING")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_READ_SPLITTING = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_READ_SPLITTING = false;
                            else
                                throw new Exception("RUN_READ_SPLITTING parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("RUN_READ_FILTERING")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_READ_FILTERING = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_READ_FILTERING = false;
                            else
                                throw new Exception("RUN_READ_FILTERING parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("FILTERING_MODE")) {
                        	try {
                        		FILTERING_MODE = FilteringMode.valueOf(valueOfParameter.trim().toUpperCase());
                        	} catch (IllegalArgumentException e) {
                        		throw new Exception("FILTERING_MODE must one of: "+Arrays.asList(FilteringMode.values())+".");
                        	}
                        } else if (nameOfParameter.equalsIgnoreCase("RUN_REFERENCE_PARTITIONING")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_REFERENCE_PARTITIONING = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_REFERENCE_PARTITIONING = false;
                            else
                                throw new Exception("RUN_REFERENCE_PARTITIONING parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("RUN_MAPPING")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_MAPPING = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_MAPPING = false;
                            else
                                throw new Exception("RUN_MAPPING parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("RUN_EXTENSION")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_EXTENSION = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_EXTENSION = false;
                            else
                                throw new Exception("RUN_EXTENSION parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("RUN_MERGE")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                RUN_MERGE = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                RUN_MERGE = false;
                            else
                                throw new Exception("RUN_MERGE parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("COMPRESS_INTERMEDIATE_FILES")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                COMPRESS_INTERMEDIATE_FILES = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                COMPRESS_INTERMEDIATE_FILES = false;
                            else
                                throw new Exception("COMPRESS_INTERMEDIATE_FILES parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("DELETE_INTERMEDIATE_FILES")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                DELETE_INTERMEDIATE_FILES = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                DELETE_INTERMEDIATE_FILES = false;
                            else
                                throw new Exception("DELETE_INTERMEDIATE_FILES parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("LENGTH_OF_READS")) {

                            try {
                                LENGTH_OF_READS = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("LENGTH_OF_READS must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MASK_OF_READS")) {

                            MASK_OF_READS = valueOfParameter;

                        } else if (nameOfParameter.equalsIgnoreCase("MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ")) {

                            try {
                                MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ")) {

                            try {
                                MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("LENGTH_OF_FIRST_PART_OF_READ")) {

                            try {
                                LENGTH_OF_FIRST_PART_OF_READ = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("LENGTH_OF_FIRST_PART_OF_READ must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("LENGTH_OF_LAST_PART_OF_READ")) {

                            try {
                                LENGTH_OF_LAST_PART_OF_READ = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("LENGTH_OF_LAST_PART_OF_READ must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING")) {

                            try {
                                MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING")) {

                            try {
                                MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION = false;
                            else
                                throw new Exception("VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION")) {

                            if (valueOfParameter.equalsIgnoreCase("YES") || valueOfParameter.equalsIgnoreCase("TRUE"))
                                MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION = true;
                            else if (valueOfParameter.equalsIgnoreCase("NO") || valueOfParameter.equalsIgnoreCase("FALSE"))
                                MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION = false;
                            else
                                throw new Exception("MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION parameter must be YES or NO");

                        } else if (nameOfParameter.equalsIgnoreCase("MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING")) {

                            try {
                                MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT")) {

                            try {
                                MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT = Integer.parseInt(valueOfParameter);
                            } catch (NumberFormatException nfe) {
                                System.out.println("MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS")) {

                            try {
                                MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS = Math.max(0, Integer.parseInt(valueOfParameter));
                            } catch (NumberFormatException nfe) {
                                System.out.println("MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS must be an integer.");
                                throw nfe;
                            }

                        } else if (nameOfParameter.equalsIgnoreCase("FILE_FILTER_SEQUENCES_FASTA")) {

                            FILE_FILTER_SEQUENCES_FASTA = new File(valueOfParameter);
                            if (!FILE_FILTER_SEQUENCES_FASTA.exists())
                                throw new Exception("FILE_FILTER_SEQUENCES_FASTA " + FILE_FILTER_SEQUENCES_FASTA.getPath() + " could not be found.");

                        } else if (nameOfParameter.equalsIgnoreCase("FILE_REFERENCE_FASTA")) {

                            FILE_REFERENCE_FASTA = new File(valueOfParameter);
                            if (!FILE_REFERENCE_FASTA.exists())
                                throw new Exception("FILE_REFERENCE_FASTA " + FILE_REFERENCE_FASTA.getPath() + " could not be found.");

                        } else if (nameOfParameter.equalsIgnoreCase("FILE_FULL_LENGTH_READS_CSFASTA")) {

                            FILE_FULL_LENGTH_READS_CSFASTA = new File(valueOfParameter);
                            if (!FILE_FULL_LENGTH_READS_CSFASTA.exists())
                                throw new Exception("FILE_FULL_LENGTH_READS_CSFASTA " + FILE_FULL_LENGTH_READS_CSFASTA.getPath() + " could not be found.");

                        } else if (nameOfParameter.equalsIgnoreCase("FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES")) {

                            FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES = new File(valueOfParameter);
                            //if (!FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES.exists())
                                //throw new Exception("FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES " + FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES.getPath() + " could not be found.");

                        } else if (nameOfParameter.equalsIgnoreCase("FOLDER_FOR_OUTPUT_FILES")) {

                            FOLDER_FOR_OUTPUT_FILES = new File(valueOfParameter);
                            if (!FOLDER_FOR_OUTPUT_FILES.exists()) {
                                System.out.println("FOLDER_FOR_OUTPUT_FILES " + FOLDER_FOR_OUTPUT_FILES.getPath() + " does not exist and will be created.");
                                FOLDER_FOR_OUTPUT_FILES.mkdirs();
                            }
                        } else if (nameOfParameter.equalsIgnoreCase("MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR")) {
                        	try { 
                        		MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR = Float.parseFloat(valueOfParameter);
                        	} catch (NumberFormatException e) {
                        		System.out.println("MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR must be an float.");
                        		throw e;
                        	}
                        } else if (nameOfParameter.equalsIgnoreCase("SCHEDULER_RESOURCE_REQUIREMENTS")) {
                        	SCHEDULER_RESOURCE_REQUIREMENTS = valueOfParameter;
                        } else if (nameOfParameter.equalsIgnoreCase("ADDITIONAL_SCHEDULER_OPTIONS")) {
                        	ADDITIONAL_SCHEDULER_OPTIONS = valueOfParameter;
                        } else if (nameOfParameter.equalsIgnoreCase("EXPERIMENTAL_MODE")) {
                        	if (Utilities.isTrue(valueOfParameter))
                        		System.setProperty(EXPERIMENTAL_MODE_SYSTEM_PROPERTY, Boolean.TRUE.toString());
                        	else if (Utilities.isFalse(valueOfParameter))
                        		System.setProperty(EXPERIMENTAL_MODE_SYSTEM_PROPERTY, Boolean.FALSE.toString());
                            else
                                throw new Exception("EXPERIMENTAL_MODE parameter must be YES or NO");
                        } else {

                            throw new Exception("Unrecognized parameter: " + nameOfParameter);
                        }


                    } else {
                        throw new Exception("The following row contains too few tab-separated columns: " + line);
                    }
                }
            }
            readerConfigurationFile.close();
	
            MASK_OF_READS = new ReadMask(MASK_OF_READS).asBitString(LENGTH_OF_READS);
            
            // check that all the required parameters were present in the config file
            if (MAX_MEMORY_PER_JOB_IN_BYTES == null)
                throw new Exception("MAX_MEMORY_PER_JOB_IN_BYTES must be specified in the config file " + fileConfiguration.getPath());
            else if (LENGTH_OF_READS == null)
                throw new Exception("LENGTH_OF_READS must be specified in the config file " + fileConfiguration.getPath());
            else if (LENGTH_OF_FIRST_PART_OF_READ == null)
                throw new Exception("LENGTH_OF_FIRST_PART_OF_READ must be specified in the config file " + fileConfiguration.getPath());
            else if (LENGTH_OF_LAST_PART_OF_READ == null)
                throw new Exception("LENGTH_OF_LAST_PART_OF_READ must be specified in the config file " + fileConfiguration.getPath());
            else if (MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING == null)
                throw new Exception("MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING must be specified in the config file " + fileConfiguration.getPath());
            else if (MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING == null)
                throw new Exception("MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING must be specified in the config file " + fileConfiguration.getPath());
            else if (MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT == null)
                throw new Exception("MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT must be specified in the config file " + fileConfiguration.getPath());
            else if (MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS == null)
                throw new Exception("MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS must be specified in the config file " + fileConfiguration.getPath());
            else if (SCORE_OF_MATCH_FOR_EXTENSION == null)
                throw new Exception("SCORE_OF_MATCH_FOR_EXTENSION must be specified in the config file " + fileConfiguration.getPath());
            else if (SCORE_OF_MISMATCH_FOR_EXTENSION == null)
                throw new Exception("SCORE_OF_MISMATCH_FOR_EXTENSION must be specified in the config file " + fileConfiguration.getPath());
			else if (RUN_READ_FILTERING && FILE_FILTER_SEQUENCES_FASTA == null)
                throw new Exception("FILE_FILTER_SEQUENCES_FASTA must be specified in the config file " + fileConfiguration.getPath());
            else if (FILE_REFERENCE_FASTA == null)
                throw new Exception("FILE_REFERENCE_FASTA must be specified in the config file " + fileConfiguration.getPath());
            else if (FILE_FULL_LENGTH_READS_CSFASTA == null)
                throw new Exception("FILE_FULL_LENGTH_READS_CSFASTA must be specified in the config file " + fileConfiguration.getPath());
            else if (FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES == null)
                throw new Exception("FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES must be specified in the config file " + fileConfiguration.getPath());
            else if (FOLDER_FOR_OUTPUT_FILES == null)
                throw new Exception("FOLDER_FOR_OUTPUT_FILES must be specified in the config file " + fileConfiguration.getPath());
            
            //Logical validation of the parameters.
            
            //Verify that the required schemas are available.
            SortedSet<File> requiredSchemaFiles = new TreeSet<File>(); 
            String leftMask = MASK_OF_READS.substring(0, LENGTH_OF_FIRST_PART_OF_READ);
            requiredSchemaFiles.add(
            	ReadMapper.getSchemaFile(leftMask, MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING, VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION, FOLDER_CONTAINING_SCHEMA_FILES)
            );
            if (LENGTH_OF_READS >= 35 && LENGTH_OF_LAST_PART_OF_READ > 0) {
            	String rightMask = "0" + MASK_OF_READS.substring(LENGTH_OF_READS - LENGTH_OF_LAST_PART_OF_READ);
            	requiredSchemaFiles.add(
            		ReadMapper.getSchemaFile(rightMask, MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING, VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION, FOLDER_CONTAINING_SCHEMA_FILES)
            	);
            };
            for ( File schemaFile : requiredSchemaFiles )
            	if (! schemaFile.exists()) throw new Exception("Schema not available: " + schemaFile.getPath());

            
            if (MAX_MEMORY_PER_JOB_IN_BYTES > 1e12)
            	warn("MAX_MEMORY_PER_JOB_IN_BYTES should be less than 1e12");
            
            if (MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR < 1 || MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR > 10)
            	throw new Exception("MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR must be between 1 and 10");
            
            if (LENGTH_OF_READS < 15) throw new Exception("LENGTH OF READS Cannot be less than 15");
            if (LENGTH_OF_READS > 100) warn("LENGTH_OF_READS is too large.");
            
            if (MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ < 1 || MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ > 100)
            	throw new Exception("MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ must be between 1 and 100");
            
            if (LENGTH_OF_FIRST_PART_OF_READ < 0)
            	throw new Exception("LENGTH_OF_FIRST_PART_OF_READ cannot be negative.");
            else if (LENGTH_OF_FIRST_PART_OF_READ > 0 && LENGTH_OF_FIRST_PART_OF_READ < 15)
            	throw new Exception("LENGTH_OF_FIRST_PART_OF_READ cannot be less than 15.");
            else if (LENGTH_OF_FIRST_PART_OF_READ > LENGTH_OF_READS)
            	throw new Exception("LENGTH_OF_FIRST_PART_OF_READ cannot be longer than LENGTH_OF_READS");
            
            if (LENGTH_OF_LAST_PART_OF_READ < 0)
            	throw new Exception("LENGTH_OF_LAST_PART_OF_READ cannot be negative.");
            else if (LENGTH_OF_LAST_PART_OF_READ > 0 && LENGTH_OF_LAST_PART_OF_READ < 15)
            	throw new Exception("LENGTH_OF_LAST_PART_OF_READ cannot be less than 15.");
            else if (LENGTH_OF_LAST_PART_OF_READ > LENGTH_OF_READS - 1)
            	throw new Exception("LENGTH_OF_LAST_PART_OF_READ cannot be longer than LENGTH_OF_READS - 1.");
            
            if (MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING < 0)
            	throw new Exception("MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING cannot be less than 0.");
            else if (MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING < 3)
            	;//OK
            else if (MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING == 3) 
            	warn("Expect long running times with MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING = 3.");
            else 
            	warn("Expect extremely long running times with MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING > 3");
            
            if (MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING < 0 )
            	throw new Exception("MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING cannot be less than 0.");
            else if (MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING < 3)
            	;//OK
            else if (MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING == 3)
            	warn("Expect long running times with MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING = 3.");
            else
            	warn("Expect extremely long running times with MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING > 3.");
            
            if (MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING < 0 )
            	throw new Exception("MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING cannot be less than 0.");
            else if (MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING < 3)
            	;//OK
            else if (MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING == 3)
            	warn("Expect long running times with MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING = 3.");
            else
            	warn("Expect extremely long running times with MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING > 3.");
            
            if (MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT < 0)
            	throw new Exception("MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT cannot be negative.");
            
            //Some config validations need the length of the reference and
            //length of the longest sequence in the reference.
            long refDbLength = 0;
            Map.Entry<String, Long> longestReferenceSequenceEntry = null;
            info("Checking reference...\n");
            FastaDatabase refDb = new FastaDatabase(FILE_REFERENCE_FASTA);
            SortedMap<String, Long> seqId2Length = refDb.getSequenceLengths();
            if (seqId2Length.size() < 1) throw new Exception("The file specified by FILE_REFERENCE_FASTA doesn't contain any sequences.");
            //Find the longest sequence.
            for (Iterator<Map.Entry<String, Long>> it = seqId2Length.entrySet().iterator(); it.hasNext();  ) {
            	Map.Entry<String, Long> entry = it.next();
            	if (longestReferenceSequenceEntry == null) longestReferenceSequenceEntry = entry;
            	else if (entry.getValue() > longestReferenceSequenceEntry.getValue()) longestReferenceSequenceEntry = entry;
            	refDbLength += entry.getValue();
            }
            
            long minimumReadToGenomeSizeRatio = 22 / 3000000000l; //Brian's suggestion.
            long readToGenomeSizeRatio = LENGTH_OF_READS / refDbLength;
            if ( readToGenomeSizeRatio < minimumReadToGenomeSizeRatio ) {
            	warn("LENGTH_OF_READS: "+ LENGTH_OF_READS + " too short for genome size: " + refDbLength + ".");
            }
            
            long mapreadsMinimumMemoryRequirement = ReadMapper.getMinimumMemoryRequirementForMapreads(longestReferenceSequenceEntry.getValue());
            if (MAX_MEMORY_PER_JOB_IN_BYTES < mapreadsMinimumMemoryRequirement)
            	warn("MAX_MEMORY_PER_JOB_IN_BYTES is set below the minimum required by mapreads.  A minimum of "+mapreadsMinimumMemoryRequirement+" is required for sequence: "+longestReferenceSequenceEntry.getKey()+", length="+longestReferenceSequenceEntry.getValue()+".");

            // End of logical validations.
            
            
            if (SCHEDULER_RESOURCE_REQUIREMENTS == null ) {
            	if (SCHEDULING_ENVIRONMENT.equalsIgnoreCase("pbs")) {
            		SCHEDULER_RESOURCE_REQUIREMENTS = JobSubmissionParameters.DEFAULT_PBS_RESOURCE_STRING;
            	} else if (SCHEDULING_ENVIRONMENT.equalsIgnoreCase("lsf")){
            		SCHEDULER_RESOURCE_REQUIREMENTS = JobSubmissionParameters.DEFAULT_LSF_RESOURCE_STRING;
            	} else if (SCHEDULING_ENVIRONMENT.equals("sge")){
            		SCHEDULER_RESOURCE_REQUIREMENTS = JobSubmissionParameters.DEFAULT_SGE_RESOURCE_STRING;
            	} else {
            		SCHEDULER_RESOURCE_REQUIREMENTS = "";
            	}
            } 
            

            System.out.println("\nVERSION\t" + VERSION);

            System.out.println("SCHEDULING_ENVIRONMENT\t" + SCHEDULING_ENVIRONMENT);
            System.out.println("NAME_OF_QUEUE\t" + NAME_OF_QUEUE);
            System.out.println("MAX_MEMORY_PER_JOB_IN_BYTES\t" + MAX_MEMORY_PER_JOB_IN_BYTES);
            System.out.println("MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES\t" + MAX_REFERENCE_SEQUENCE_PER_JOB_IN_BASES);
            System.out.println("NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE\t" + NUMBER_OF_READS_PER_MERGE_JOB_OVERRIDE);
            System.out.println("NUMBER_OF_READS_FILE_PARTITIONS\t" + NUMBER_OF_READS_FILE_PARTITIONS);
            System.out.println("JOBS_RERUNNABLE\t" + JOBS_RERUNNABLE);

            System.out.println("RUN_READ_SPLITTING\t" + RUN_READ_SPLITTING);
            System.out.println("RUN_READ_FILTERING\t" + RUN_READ_FILTERING);
            System.out.println("RUN_REFERENCE_PARTITIONING\t" + RUN_REFERENCE_PARTITIONING);
            System.out.println("RUN_MAPPING\t" + RUN_MAPPING);
            System.out.println("RUN_EXTENSION\t" + RUN_EXTENSION);
            System.out.println("RUN_MERGE\t" + RUN_MERGE);

            System.out.println("FILTERING_MODE\t" + FILTERING_MODE);

            System.out.println("COMPRESS_INTERMEDIATE_FILES\t" + COMPRESS_INTERMEDIATE_FILES);
            System.out.println("DELETE_INTERMEDIATE_FILES\t" + DELETE_INTERMEDIATE_FILES);

            System.out.println("LENGTH_OF_READS\t" + LENGTH_OF_READS);
            System.out.println("MASK_OF_READS\t" + MASK_OF_READS);
            System.out.println("MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ\t" + MIN_MAPPING_LOCATIONS_REQUIRED_BEFORE_REPORTING_READ);
            System.out.println("MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ\t" + MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ);
            System.out.println("VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION\t" + VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE_FOR_MAPPING_AND_EXTENSION);
            System.out.println("MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION\t" + MATCHES_TO_IUPACS_VALID_FOR_MAPPING_AND_EXTENSION);
            System.out.println("LENGTH_OF_FIRST_PART_OF_READ\t" + LENGTH_OF_FIRST_PART_OF_READ);
            System.out.println("LENGTH_OF_LAST_PART_OF_READ\t" + LENGTH_OF_LAST_PART_OF_READ);
            System.out.println("MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING\t" + MAX_MISMATCHES_IN_FIRST_PART_OF_READ_FOR_MAPPING);
            System.out.println("MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING\t" + MAX_MISMATCHES_IN_LAST_PART_OF_READ_FOR_MAPPING);
            System.out.println("MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING\t" + MAX_MISMATCHES_IN_READ_PARTS_FOR_READ_FILTERING);
            System.out.println("MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT\t" + MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT);
            System.out.println("MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS\t" + MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS);
            System.out.println("SCORE_OF_MATCH_FOR_EXTENSION\t" + SCORE_OF_MATCH_FOR_EXTENSION);
            System.out.println("SCORE_OF_MISMATCH_FOR_EXTENSION\t" + SCORE_OF_MISMATCH_FOR_EXTENSION);
            if (Utilities.isTrue(System.getProperty(EXPERIMENTAL_MODE_SYSTEM_PROPERTY)))
            	System.out.println("EXPERIMENTAL_MODE\ttrue\n");
            
            System.out.println("FOLDER_INSTALLATION_ROOT\t" + FOLDER_INSTALLATION_ROOT.getPath());
            System.out.println("FILE_FILTER_SEQUENCES_FASTA\t" + FILE_FILTER_SEQUENCES_FASTA.getPath());
            System.out.println("FILE_REFERENCE_FASTA\t" + FILE_REFERENCE_FASTA.getPath());
            System.out.println("FILE_FULL_LENGTH_READS_CSFASTA\t" + FILE_FULL_LENGTH_READS_CSFASTA.getPath());
            System.out.println("FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES\t" + FOLDER_FOR_TEMPORARY_FILES_ON_COMPUTE_NODES.getPath());
            System.out.println("FOLDER_FOR_OUTPUT_FILES\t" + FOLDER_FOR_OUTPUT_FILES.getPath());

            if (MASK_OF_READS != null && MASK_OF_READS.length() != LENGTH_OF_READS)
                throw new Exception("MASK_OF_READS must be the same length as LENGTH_OF_READS (" + MASK_OF_READS.length() + " != " + LENGTH_OF_READS + ").");

            if (MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ > 10)
                System.out.println("WARNING: MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ of " + MAX_MAPPING_LOCATIONS_ALLOWED_BEFORE_NOT_REPORTING_READ + " is greater than our highest tested value of 10 and could potentially lead to performance issues.");            
            if (MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT < 24)
                System.out.println("WARNING: MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT of " + MIN_ALIGNMENT_SCORE_FOR_REPORTING_ALIGNMENT + " is less than our lowest tested value of 24 and could potentially lead to performance issues.");
            if (MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS > 4)
                System.out.println("WARNING: MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS of " + MIN_GAP_IN_ALIGNMENT_SCORE_TO_SECOND_BEST_ALIGNMENT_FOR_UNIQUENESS + " is less than our lowest tested value and may be overly stringent.");
            if (LENGTH_OF_FIRST_PART_OF_READ < 25)
                System.out.println("WARNING: LENGTH_OF_FIRST_PART_OF_READ of " + LENGTH_OF_FIRST_PART_OF_READ + " is less than our lowest tested value of 25 and could potentially lead to performance issues.");
            if (LENGTH_OF_LAST_PART_OF_READ < 25 && LENGTH_OF_LAST_PART_OF_READ > 1)
                System.out.println("WARNING: LENGTH_OF_LAST_PART_OF_READ of " + LENGTH_OF_LAST_PART_OF_READ + " is less than our lowest tested value of 25 and could potentially lead to performance issues.");


        } catch (Exception e) {
            System.out.println("Error encontered while loading configuration file.");
            throw e;
        }


    }

    /**
     * Current behavior is to throw a generic exception if one or more mapping jobs failed.
     *
     * @param setOfScriptOutputFilePaths
     * @throws Exception
     */
    public static void checkThatAllJobsCompletedSuccessfully(HashSet<String> setOfScriptOutputFilePaths) throws Exception {

        ArrayList<File> logOfScriptOutputFilesFromPreviousRun = new ArrayList<File>();
        Iterator<String> iteratorOverScriptOutputFilePaths = setOfScriptOutputFilePaths.iterator();
        boolean noMissingFiles = true;
        while (iteratorOverScriptOutputFilePaths.hasNext()) {
            File fileScriptOutput = new File( iteratorOverScriptOutputFilePaths.next());
            if (!fileScriptOutput.exists()) {
                if (noMissingFiles) {
                    noMissingFiles = false;
                    System.out.println(" Fail\nThe following script output files could not be found:");
                }
                System.out.println(" -->" + fileScriptOutput.getPath());
            }
            logOfScriptOutputFilesFromPreviousRun.add(fileScriptOutput);
        }

        if (!noMissingFiles)
            throw new Exception("One or more job output files could not be found and therefore these jobs did not complete or failed.");

        ArrayList<File> listOfScriptOutputFilesForFailedJobs = ClusterInterface.getListOfScriptOutputFilesIndicatingFailure(logOfScriptOutputFilesFromPreviousRun);
        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
            System.out.println(" Fail\nThe following jobs failed:");
            Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
            while (iteratorOverFailFiles.hasNext())
                System.out.println(" -->" + ( iteratorOverFailFiles.next()).getPath());

            throw new Exception("One or more jobs failed.");
        }
    }

    @SuppressWarnings("unchecked")
    public static ArrayList<Integer>[] loadPartitionTable(File filePartitionTable) throws IOException {
        ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
        int maxRefPartitionIndex = getMaxPartitionIndexFromPartitionTableFile(filePartitionTable);
        listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = new ArrayList[maxRefPartitionIndex +1];
        
        for (int i=0; i<listsOfFullRefSequenceNumbersOrderedForEachReferencePartition.length; i++)
        	listsOfFullRefSequenceNumbersOrderedForEachReferencePartition[i] = new ArrayList<Integer>();
        
        BufferedReader readerReferencePartitionTable = null;
        try {
	        readerReferencePartitionTable = new BufferedReader(new FileReader(filePartitionTable));
	        for (String line = readerReferencePartitionTable.readLine(); line != null; line = readerReferencePartitionTable.readLine()) {
	        	if (line.trim().startsWith("#")) continue;
	            String tokens[] = line.split("\t");
	            int partitionNum = Integer.parseInt(tokens[0]);
	            int seqNum = Integer.parseInt(tokens[2]);
	            listsOfFullRefSequenceNumbersOrderedForEachReferencePartition[partitionNum].add(seqNum);
	        }
        } finally {
        	if (readerReferencePartitionTable != null) readerReferencePartitionTable.close();
        }
        return listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
    }

    private static int getMaxPartitionIndexFromPartitionTableFile(File filePartitionTable) throws IOException {
        BufferedReader readerReferencePartitionTable = null;
        try {
        	readerReferencePartitionTable = new BufferedReader(new FileReader(filePartitionTable));
	        int maxRefPartitionIndex = 0;
	        for (String line = readerReferencePartitionTable.readLine(); line != null; line=readerReferencePartitionTable.readLine()) {
	        	if (line.trim().startsWith("#")) continue; //Drop the comments.
	        	if (line.startsWith("refPartitionIndex")) continue; //Drop the header line
	            maxRefPartitionIndex = Math.max(maxRefPartitionIndex, Integer.parseInt(line.split("\t")[0]));
	        }
	        return maxRefPartitionIndex;
        } finally {
        	if (readerReferencePartitionTable != null) readerReferencePartitionTable.close();
        }
    }
    
    public static void info(Object msg) {
    	System.out.println(msg);
    }
    
    public static void warn(Object msg) {
    	System.out.println("WARNING! " + msg);
    }
    
    public static File getAppRoot() {
    	return new File(System.getProperty(WT_HOME_SYSTEM_PROPERTY));
    }
    
    public static File getMapreadsExe() {
    	return new File(getAppRoot(), "pkg/mapreads");
    }
    
    public static File getExtendMappedReadsPy() {
    	return new File(getAppRoot(), "pkg/extendMappedReads.py");
    }
    	
    public static File getSchemaDir() {
    	return new File(getAppRoot(), "etc/schemas/");
    }
    
    public static File getLibDir() {
    	return new File(getAppRoot(), "etc/lib");
    }
    	
    public static File getWholeTranscriptomeJar() {
    	return new File(getAppRoot(), "pkg/WholeTranscriptome.jar");
    }
}

