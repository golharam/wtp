package com.lifetechnologies.solid.wt;

import java.io.File;
import java.util.ArrayList;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;

/**
 * User: tuchbb
 * Date: Dec 1, 2008
 * Time: 10:50:43 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExonFinderOptimizer {

    public static void main(String args[]) throws Exception {

        int widthsOfMovingWindow[] = { 25, 50, 75, 100, 125, 150, 200 };
        int minCoveragesToCallAsPartOfExon[] = { 100, 250, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 4000, 5000 };
        double minOverlapsForCountingPredictedExonAsTruePositive[] = { 0.5 };   //0.01, 0.25, 0.5 };
        int minAlignmentScore = 24;
        double minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon[] = { 0.0, 0.1, 0.25, 0.5 };
        JobSubmissionParameters params = new JobSubmissionParameters();
        params.setEnvironment("pbs");
        params.setQueueName("highmem");
        params.setResourceString(JobSubmissionParameters.DEFAULT_PBS_RESOURCE_STRING);
        params.setMemoryRequirement(new Double(1.5e9).longValue());
        params.setRerunnable(false);
        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
        File fileListOfSortedMAX = new File("/home/tuchbb/projects/mayo/mapping_results/listOfSortedMAXFiles.txt");
        File fileReferenceFasta = new File("/share/reference/genomes/chromFa/human.fa");
        File fileExonAnnotationGFF = new File("/home/tuchbb/data/species/h_sapiens/ucsc_hg18_081115/exon_annotation.filtered.081115.gff");
        //File folderForOutput = new File("/home/tuchbb/projects/MAQC/mapping_results/HBR_polyA/match_genome_25.2_30.2_mask5/output/exonFinder/");
        File folderForOutput = new File("/data/results/mayo/results/exonFinder/");

        for (int indexWidthOfMovingWindow = 0; indexWidthOfMovingWindow < widthsOfMovingWindow.length; indexWidthOfMovingWindow++)
            for (int indexMinCoverage = 0; indexMinCoverage < minCoveragesToCallAsPartOfExon.length; indexMinCoverage++)
                for (int indexMinTrimCoverage = 0; indexMinTrimCoverage < minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon.length; indexMinTrimCoverage++)
                    for (int indexMinOverlap = 0; indexMinOverlap < minOverlapsForCountingPredictedExonAsTruePositive.length; indexMinOverlap++) {
                        String prefixFilePath = folderForOutput.getPath() + "/exonFinder." +  widthsOfMovingWindow[indexWidthOfMovingWindow] + "_" + minCoveragesToCallAsPartOfExon[indexMinCoverage]
                                                    + "_" + minOverlapsForCountingPredictedExonAsTruePositive[indexMinOverlap] + "_" + minAlignmentScore
                                                    + "_" + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon[indexMinTrimCoverage] + ".";

                        ArrayList<String> commandStrings = new ArrayList<String>();
                        commandStrings.add("java -cp /home/tuchbb/bin/java/geneExpression/bin/ -Xmx" + (int)Math.ceil(clusterInterface.getJobSubmissionParameters().getMemoryRequirement() / Constants.BYTES_PER_MEGABYTE) + "m " +
                                            ExonFinder.class.getName() + " " + widthsOfMovingWindow[indexWidthOfMovingWindow] + " " + minCoveragesToCallAsPartOfExon[indexMinCoverage] + " " +
                                            minOverlapsForCountingPredictedExonAsTruePositive[indexMinOverlap] + " " + minAlignmentScore + " " + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon[indexMinTrimCoverage] + " " +
                                            fileListOfSortedMAX.getPath() + " " + fileReferenceFasta.getPath() + " " + fileExonAnnotationGFF.getPath() + " " + folderForOutput.getPath() +
                                            " > " + prefixFilePath + "output.txt" );
                        clusterInterface.executeJob(new File(prefixFilePath + "sh"), commandStrings);

                    }
        File fileJobRemoval = new File(folderForOutput.getPath() + "/remove_jobs.sh");
        clusterInterface.writeMasterJobRemovalFileFromLog(fileJobRemoval);
        fileJobRemoval.setExecutable(true);
    }
}
