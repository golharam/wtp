package com.lifetechnologies.solid.wt;

import java.util.*;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * User: tuchbb
 * Date: Dec 1, 2008
 * Time: 1:53:40 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ContiguousGenomicFeatureAlignmentCounter {

    public static void main(String args[]) throws Exception {

        int minAlignmentScore = Integer.parseInt(args[0]);
        int numberOfBasesToLengthenFeaturesBy = Integer.parseInt(args[1]);
        File fileListOfSortedMAXFiles = new File(args[2]);
        //SortedMAXFilesBufferedReader readerSortedMAXFiles = new SortedMAXFilesBufferedReader(filesSortedMax);
        FastaDatabase fastaDatabaseOfReference = new FastaDatabase(new File(args[3]));
        File fileExonAnnotationGFF = new File(args[4]);
        File fileHistogramOfCoveragePDF = new File(args[5]);
        String suffixCountsFile = args[6];  //exonCounts.txt

        Histogram histogramOfCoverageByFeature = new Histogram(51, 2, 0.0, "coverage");
        HashSet<String> setOfMaxFilesPaths = TextFileUtilities.loadSetFromFile(fileListOfSortedMAXFiles, "\t", 0);
        Iterator<String> iteratorOverMAXFilePaths = setOfMaxFilesPaths.iterator();
        while (iteratorOverMAXFilePaths.hasNext()) {
            File filesSortedMAX[] = { new File(iteratorOverMAXFilePaths.next()) };
            File fileCountsOutput = new File(filesSortedMAX[0].getPath().replaceAll("csfasta$", fileExonAnnotationGFF.getName().replaceAll(".gff$", "") + "." + suffixCountsFile) );

            System.out.println("Processing " + filesSortedMAX[0].getPath());
            System.out.println("Output will be " + fileCountsOutput.getPath());
            count(minAlignmentScore, numberOfBasesToLengthenFeaturesBy, histogramOfCoverageByFeature, filesSortedMAX[0].getName(), false,
                    new SortedMAXFilesBufferedReader(filesSortedMAX, new ExtendedReadMappingByPositionComparator(), false), fastaDatabaseOfReference, fileExonAnnotationGFF, fileCountsOutput);
        }
        System.out.println("\n" + histogramOfCoverageByFeature);
        histogramOfCoverageByFeature.toPDF(fileHistogramOfCoveragePDF, true, false, false, null);

    }

    private static void count(int minAlignmentScore, int numberOfBasesToLengthenFeaturesBy, Histogram histogramOfCoverageByFeature, String nameOfHistogramSeries, boolean printNumberOfUniqueAlignmentStartPositions,
                              SortedMAXFilesBufferedReader readerSortedMAXFiles, FastaDatabase fastaDatabaseOfReference,
                              File fileFeatureAnnotationGFF, File fileCountsOutput) throws Exception, CloneNotSupportedException {

        List<String> listOfHeaders = fastaDatabaseOfReference.getListOfHeaders(true, false);

        ContiguousGenomicFeatureList listOfContiguousGenomicFeatures = new ContiguousGenomicFeatureList(fileFeatureAnnotationGFF, 1000, numberOfBasesToLengthenFeaturesBy, true, true, false, true);

        TreeMap<ContiguousGenomicFeature, HashSet<Integer>> mapFeatureToAlignmentStartPositions = new TreeMap<ContiguousGenomicFeature, HashSet<Integer>>(new ContiguousGenomicFeatureComparator());
        Iterator<ContiguousGenomicFeature> iteratorOverAllFeatures = listOfContiguousGenomicFeatures.getSetOfAllFeatures().iterator();
        while (iteratorOverAllFeatures.hasNext()) {
            ContiguousGenomicFeature contiguousGenomicFeature = iteratorOverAllFeatures.next();
            contiguousGenomicFeature.setValueA(0.0);
            contiguousGenomicFeature.setValueB(0.0);
            mapFeatureToAlignmentStartPositions.put(contiguousGenomicFeature, null);
        }

        int countReads = 0;
        int countReadsWithMinAlignScore = 0;
        int countReadsWithMinAlignScoreInFeature = 0;

        String headerOfCurrentReferenceSequence;
        String headerOfLastReferenceSequence = "";

        ExtendedReadMapping extendedReadMappingCurrent;
        while ((extendedReadMappingCurrent = readerSortedMAXFiles.nextEntry()) != null) {

            headerOfCurrentReferenceSequence = listOfHeaders.get(extendedReadMappingCurrent.getIndexOfMatchingReferenceSequence() -1);
            if (!headerOfCurrentReferenceSequence.equalsIgnoreCase(headerOfLastReferenceSequence))
                System.out.println("Processing reference sequence:\t" + headerOfCurrentReferenceSequence);
            headerOfLastReferenceSequence = headerOfCurrentReferenceSequence;

            countReads++;

            if (extendedReadMappingCurrent.getScore() >= minAlignmentScore) {

                countReadsWithMinAlignScore++;

                int positionAbsoluteWhereAlignmentStartsInReferenceOneBased = Math.abs(extendedReadMappingCurrent.getPositionOfAlignmentStartInReferenceSequence()) +1;
                int positionAbsoluteWhereAlignmentEndsInReferenceOneBased = positionAbsoluteWhereAlignmentStartsInReferenceOneBased + extendedReadMappingCurrent.getLengthOfAlignment();

                Strand strand = (extendedReadMappingCurrent.getPositionOfAlignmentStartInReferenceSequence() >= 0) ? Strand.POSITIVE : Strand.NEGATIVE;
                HashSet<ContiguousGenomicFeature> features = listOfContiguousGenomicFeatures.getSetOfFeaturesOverlapping(headerOfCurrentReferenceSequence, strand, positionAbsoluteWhereAlignmentStartsInReferenceOneBased, positionAbsoluteWhereAlignmentEndsInReferenceOneBased, true).getFeaturesThatOverlapThisOne();

                Iterator<ContiguousGenomicFeature> iteratorOverFeatures = features.iterator();
                while (iteratorOverFeatures.hasNext()) {
                    ContiguousGenomicFeature feature = iteratorOverFeatures.next();
                    HashSet<Integer> setOfPositionsWhereAlignmentsStart = mapFeatureToAlignmentStartPositions.get(feature);
                    if (setOfPositionsWhereAlignmentsStart == null) {
                        setOfPositionsWhereAlignmentsStart = new HashSet<Integer>();
                        mapFeatureToAlignmentStartPositions.put(feature, setOfPositionsWhereAlignmentsStart);
                    }
                    feature.setValueA(feature.getValueA() + 1.0);
                    feature.setValueB(feature.getValueB() + extendedReadMappingCurrent.getLengthOfAlignment());
                    setOfPositionsWhereAlignmentsStart.add(positionAbsoluteWhereAlignmentStartsInReferenceOneBased);

                }

                if (features.size() > 0)
                    countReadsWithMinAlignScoreInFeature++;
            }
        }


        BufferedWriter writerCountsOutputFile = new BufferedWriter(new FileWriter(fileCountsOutput));
        writerCountsOutputFile.write("#countReads\t" + countReads);
        writerCountsOutputFile.newLine();
        writerCountsOutputFile.write("#countReadsWithMinAlignScore\t" + countReadsWithMinAlignScore);
        writerCountsOutputFile.newLine();
        writerCountsOutputFile.write("#countReadsWithMinAlignScoreInFeature\t" + countReadsWithMinAlignScoreInFeature);
        writerCountsOutputFile.newLine();

        int countFeaturesWithAtLeastOneAlignment = 0;
        int countFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions = 0;
        int countFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions = 0;
        double totalCoverageOfFeatures = 0;
        double totalLengthOfFeatures = 0;
        Iterator<ContiguousGenomicFeature> iteratorOverFeatures = mapFeatureToAlignmentStartPositions.keySet().iterator();
        while (iteratorOverFeatures.hasNext()) {
            ContiguousGenomicFeature feature =  iteratorOverFeatures.next();
            HashSet<Integer> setOfPositionsWithAlignmentStart = mapFeatureToAlignmentStartPositions.get(feature);
            if (setOfPositionsWithAlignmentStart != null) {
                countFeaturesWithAtLeastOneAlignment++;
                if (setOfPositionsWithAlignmentStart.size() > 1)
                    countFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions++;
                if (setOfPositionsWithAlignmentStart.size() > 3)
                    countFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions++;
            }
            totalCoverageOfFeatures += feature.getValueB();
            totalLengthOfFeatures += feature.getLength();
            histogramOfCoverageByFeature.addObservationToHistogramSeries(nameOfHistogramSeries, feature.getValueB() / feature.getLength());

            if (printNumberOfUniqueAlignmentStartPositions) {
                if (setOfPositionsWithAlignmentStart == null)   feature.setValueB(0.0);
                else    feature.setValueB((double)setOfPositionsWithAlignmentStart.size());
            } else
                feature.setValueB(null);

            writerCountsOutputFile.write(feature.toString());
            writerCountsOutputFile.newLine();
        }
        writerCountsOutputFile.close();

        double averageCoverageOfFeatures = totalCoverageOfFeatures / totalLengthOfFeatures;

        System.out.println("countReads:\t" + countReads);
        System.out.println("countReadsWithMinAlignScore:\t" + countReadsWithMinAlignScore);
        System.out.println("countReadsWithMinAlignScoreInFeature:\t" + countReadsWithMinAlignScoreInFeature);
        System.out.println("countFeatures:\t" + mapFeatureToAlignmentStartPositions.size());
        System.out.println("countFeaturesWithAtLeastOneAlignment:\t" + countFeaturesWithAtLeastOneAlignment);
        System.out.println("countFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions:\t" + countFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions);
        System.out.println("countFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions:\t" + countFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions);
        System.out.println("averageCoverageOfFeatures:\t" + averageCoverageOfFeatures);

        System.out.println("\ncountReads\tcountReadsWithMinAlignScore\tcountReadsWithMinAlignScoreInFeature\tcountFeatures\tcountFeaturesWithAtLeastOneAlignment\tcountFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions\tcountFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions\taverageCoverageOfFeatures");
        System.out.println(countReads + "\t" + countReadsWithMinAlignScore + "\t" + countReadsWithMinAlignScoreInFeature + "\t" + mapFeatureToAlignmentStartPositions.size() + "\t"
                            + countFeaturesWithAtLeastOneAlignment  + "\t" + countFeaturesWithAtLeastTwoAlignmentsWithUniqueStartPositions + "\t" + countFeaturesWithAtLeastFourAlignmentsWithUniqueStartPositions + "\t" + averageCoverageOfFeatures);
        System.out.println("\n");



    }
}
