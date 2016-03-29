package com.lifetechnologies.solid.wt;

import java.util.TreeMap;
import java.util.Iterator;
import java.util.ArrayList;
import java.io.*;
import java.text.DecimalFormat;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;

/**
 * User: tuchbb
 * Date: Nov 26, 2008
 * Time: 9:58:23 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class MAXFileSorter {

    private static final String PATH_TO_JAVA_EXE = "java";

    public static void main(String args[]) throws Exception {

        int minAlignmentScoreForReportingAlignment = Integer.parseInt(args[0]);
        int minAlignmentsWithMinScoreRequiredBeforeReportingRead = Integer.parseInt(args[1]);
        int maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead = Integer.parseInt(args[2]);
        int minScoreGapToSecondBestAlignment = Integer.parseInt(args[3]);

        File folderInstallationRoot = new File(args[4]);    //"/home/tuchbb/bin/java/geneExpression/");
        File fileMAXUnsorted = new File(args[5]);   //"/home/tuchbb/projects/MAQC/mapping_results/HBR_polyA/match_genome_25.2_30.2_mask5/output/73_20080930_1_Ambion_HBR_HBR_F3.max.merged.filtered.csfasta");
        File folderTemp = new File(args[6]);    //"/home/tuchbb/projects/MAQC/mapping_results/HBR_polyA/match_genome_25.2_30.2_mask5/output/sort_tmp/");
        if (!folderTemp.exists())
            folderTemp.mkdirs();
        File fileMAXSortedOutput = new File(args[7]);   //"/home/tuchbb/projects/MAQC/mapping_results/HBR_polyA/match_genome_25.2_30.2_mask5/output/73_20080930_1_Ambion_HBR_HBR_F3.max.merged.filtered.sorted_24_1_10_0.csfasta");

        String nameOfClusterQueue = args[8];
	    long numberOfBytesOfRAMAllocatedPerJob = Long.parseLong(args[9]);
        String schedulingEnvironment = args.length > 10 ? args[10] : "pbs";

        String stringClassPath = folderInstallationRoot.getPath() + "/bin/" + ":" + folderInstallationRoot.getPath() + "/lib/" + "/*";

        JobSubmissionParameters params = new JobSubmissionParameters();
        params.setEnvironment(schedulingEnvironment);
        params.setQueueName(nameOfClusterQueue);
        params.setRerunnable(false);
        params.setMemoryRequirement(numberOfBytesOfRAMAllocatedPerJob);
        if (schedulingEnvironment.equals("pbs")) {
        	params.setResourceString(JobSubmissionParameters.DEFAULT_PBS_RESOURCE_STRING);
        } else if (schedulingEnvironment.equals("lsf")){
        	params.setResourceString(JobSubmissionParameters.DEFAULT_LSF_RESOURCE_STRING);
        } else {
        	throw new Exception("Unknown scheduling environment: "+schedulingEnvironment);
        }
        
        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);

        sort(minAlignmentScoreForReportingAlignment,
                minAlignmentsWithMinScoreRequiredBeforeReportingRead,
                maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead,
                minScoreGapToSecondBestAlignment,
                clusterInterface, stringClassPath, fileMAXUnsorted, folderTemp, fileMAXSortedOutput);

    }


    public static void sort(int minAlignmentScoreForReportingAlignment,
                            int minAlignmentsWithMinScoreRequiredBeforeReportingRead,
                            int maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead,
                            int minScoreGapToSecondBestAlignment,
                            ClusterInterface clusterInterface,
                            String stringClassPathForJava,
                            File fileMAXUnsorted,
                            File folderTemp,
                            File fileMAXSortedOutput) throws Exception {


        // index the MAX file        
        int countOfReads = 0;
        long countOfCharactersRead = 0;
        BufferedReader readerMAXFile = new BufferedReader(new FileReader(fileMAXUnsorted));
        File fileIndexOfMAXFile = new File(fileMAXUnsorted.getPath() + ".index");
        BufferedWriter writerIndexOfMAXFile = new BufferedWriter(new FileWriter(fileIndexOfMAXFile));
        String line;
        while ((line = readerMAXFile.readLine()) != null) {
            if (countOfReads % 1000000 == 0) {
                writerIndexOfMAXFile.write(countOfReads + "\t" + new DecimalFormat("0").format(countOfCharactersRead));
                writerIndexOfMAXFile.newLine();
                writerIndexOfMAXFile.flush();
            }
            countOfCharactersRead += line.length() +1;
            countOfCharactersRead += readerMAXFile.readLine().length() +1;
            countOfReads++;
        }
        readerMAXFile.close();
        writerIndexOfMAXFile.close();        

        // sort partitions
        int numberOfReadsPerPartition = (int) (clusterInterface.getJobSubmissionParameters().getMemoryRequirement() / 4000.0);
        int numberOfPartitions = (int)Math.ceil(countOfReads / numberOfReadsPerPartition);
        File filesSortedOutputForEachReadsFilePartition[] = new File[numberOfPartitions];
        for (int indexReadsFilePartition = 0; indexReadsFilePartition < filesSortedOutputForEachReadsFilePartition.length; indexReadsFilePartition++) {

            ArrayList<String> commandStrings = new ArrayList<String>();

            int indexOfFirstReadInPartition = numberOfReadsPerPartition * indexReadsFilePartition;
            int indexOfLastReadInPartition = Math.min(numberOfReadsPerPartition * (indexReadsFilePartition +1) -1, countOfReads -1);
            filesSortedOutputForEachReadsFilePartition[indexReadsFilePartition] = new File(folderTemp.getPath() + "/reads_" + indexOfFirstReadInPartition + "_" + indexOfLastReadInPartition + ".sorted.csfasta");

            String command = PATH_TO_JAVA_EXE + " -classpath " + stringClassPathForJava + " -Xmx" + (int)Math.ceil(clusterInterface.getJobSubmissionParameters().getMemoryRequirement()/ Constants.BYTES_PER_MEGABYTE) + "m" +
                                " " + MAXFileSorterSingleNode.class.getName() +
                                " " + indexOfFirstReadInPartition +
                                " " + indexOfLastReadInPartition +
                                " " + minAlignmentScoreForReportingAlignment +
                                " " + minAlignmentsWithMinScoreRequiredBeforeReportingRead +
                                " " + maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead +
                                " " + minScoreGapToSecondBestAlignment +
                                " " + fileMAXUnsorted.getPath() +
                                " " + fileIndexOfMAXFile.getPath() +
                                " " + filesSortedOutputForEachReadsFilePartition[indexReadsFilePartition].getPath();

            commandStrings.add(command);

            File fileScript = new File(filesSortedOutputForEachReadsFilePartition[indexReadsFilePartition].getPath().replaceAll(".csfasta$", ".sh"));

            clusterInterface.executeJob(fileScript, commandStrings);
        }

        File fileJobRemoval = new File(folderTemp.getPath() + "/remove_jobs.sh");
        clusterInterface.writeMasterJobRemovalFileFromLog(fileJobRemoval);
        fileJobRemoval.setExecutable(true);

        System.out.print("\nWaiting for sort jobs to finish...");
        while (!clusterInterface.checkIfLoggedJobsComplete())
            Thread.sleep(5000);

        System.out.println("OK");

        ArrayList<File> listOfScriptOutputFilesForFailedJobs = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
        if (listOfScriptOutputFilesForFailedJobs.size() > 0) {
            System.out.println("The following jobs failed:");
            Iterator<File> iteratorOverFailFiles = listOfScriptOutputFilesForFailedJobs.iterator();
            while (iteratorOverFailFiles.hasNext())
                System.out.println((iteratorOverFailFiles.next()).getPath());

            throw new Exception("One or more sort jobs failed.");
        }


        System.out.println("\n -->Merging sorted MAX files, by reads file partition.");
        BufferedWriter writer = new BufferedWriter(new FileWriter(fileMAXSortedOutput));
        BufferedReader readersSortedOutputFilesByReadsFilePartition[] = new BufferedReader[filesSortedOutputForEachReadsFilePartition.length];
        for (int i = 0; i < filesSortedOutputForEachReadsFilePartition.length; i++)
            readersSortedOutputFilesByReadsFilePartition[i] = new BufferedReader(new FileReader(filesSortedOutputForEachReadsFilePartition[i]));

        TreeMap<ExtendedReadMapping, Integer> mapExtendedReadMappingToPartitionIndex = new TreeMap<ExtendedReadMapping, Integer>(new ExtendedReadMappingByPositionComparator());
        String lines[] = new String[filesSortedOutputForEachReadsFilePartition.length];
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

    }



}
