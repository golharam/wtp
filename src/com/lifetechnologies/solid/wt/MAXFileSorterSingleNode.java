package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.TreeMap;
import java.util.Iterator;

/**
 * User: tuchbb
 * Date: Nov 26, 2008
 * Time: 3:42:47 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class MAXFileSorterSingleNode {

    public static void main(String args[]) throws IOException {

        int indexOfFirstReadInPartition = Integer.parseInt(args[0]);
        int indexOfLastReadInPartition = Integer.parseInt(args[1]);
        int minAlignmentScoreForReportingAlignment = Integer.parseInt(args[2]);
        int minAlignmentsWithMinScoreRequiredBeforeReportingRead = Integer.parseInt(args[3]);
        int maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead = Integer.parseInt(args[4]);
        int minScoreGapToSecondBestAlignment = Integer.parseInt(args[5]);

        File fileMAXUnsorted = new File(args[6]);
        File fileIndexOfMAXUnsorted = new File(args[7]);
        File fileMAXSortedOutput = new File(args[8]);

        sort(indexOfFirstReadInPartition,
                indexOfLastReadInPartition,
                minAlignmentScoreForReportingAlignment,
                minAlignmentsWithMinScoreRequiredBeforeReportingRead,
                maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead,
                minScoreGapToSecondBestAlignment,
                fileMAXUnsorted,
                fileIndexOfMAXUnsorted,
                fileMAXSortedOutput);

    }        

    private static void sort(int indexOfFirstReadToSort,
                                int indexOfLastReadToSort,
                                int minAlignmentScoreForReportingAlignment,
                                int minAlignmentsWithMinScoreRequiredBeforeReportingRead,
                                int maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead,
                                int minScoreGapToSecondBestAlignment,
                                File fileMAXUnsorted,
                                File fileIndexOfMAXUnsorted,
                                File fileMAXSortedOutput) throws IOException {

        TreeMap<ExtendedReadMapping, Read> mapSortedReadMappingToMergedMaxFileEntry = new TreeMap<ExtendedReadMapping, Read>(new ExtendedReadMappingByPositionComparator());

        int countOfReadsProcessed = 0;
        int countOfReadsWithNumberOfAlignmentsInRequiredRange = 0;

        BufferedRandomAccessFile bufferedRandomAccessFileMAX = new BufferedRandomAccessFile(fileMAXUnsorted, "r");

        // advance to the (indexOfFirstReadToSort +1)th read
        int indexOfCurrentRead = 0;
        TreeMap<String, Long> mapReadIndexToCharacterIndex = TextFileUtilities.loadStringToLongMapFromFile(fileIndexOfMAXUnsorted, 0, 1, "\t", 0);

        int indexOfReadToSeekTo = indexOfFirstReadToSort;
        long indexOfCharToSeekTo = -1;
        while (indexOfCharToSeekTo < 0 && indexOfReadToSeekTo >= 0) {
            if (mapReadIndexToCharacterIndex.containsKey("" + indexOfReadToSeekTo))
                indexOfCharToSeekTo = mapReadIndexToCharacterIndex.get("" + indexOfReadToSeekTo);
            else
                indexOfReadToSeekTo--;
        }

        int numberOfCharsPerSkip = 2000000000;
        int numberOfSkips = (int)Math.ceil(indexOfCharToSeekTo / numberOfCharsPerSkip);
        for (int i = 0; i < numberOfSkips -1; i++)
            bufferedRandomAccessFileMAX.skipBytes(numberOfCharsPerSkip);

        bufferedRandomAccessFileMAX.skipBytes((int)(indexOfCharToSeekTo % numberOfCharsPerSkip));

        indexOfCurrentRead = indexOfReadToSeekTo;
        String line = bufferedRandomAccessFileMAX.readNextLine(); // header

        System.out.println("Skipped " + indexOfCurrentRead + " reads by seeking past " + indexOfCharToSeekTo + " characters.");
        System.out.println(" Current read is:\t" + line);


        while (indexOfCurrentRead++ < indexOfFirstReadToSort) {
            bufferedRandomAccessFileMAX.readNextLine();    // sequence
            line = bufferedRandomAccessFileMAX.readNextLine(); // header
            if (indexOfCurrentRead % 10000000 == 0)
                System.out.println("Skipped " + indexOfCurrentRead + " reads.");
        }

        System.out.println(" Current read is:\t" + line);


        // each loop iteration deals with one read
        while (line != null && countOfReadsProcessed < (indexOfLastReadToSort - indexOfFirstReadToSort +1)) {

            countOfReadsProcessed++;

            if (countOfReadsProcessed % 1000000 == 0)
                System.out.println("Processed " + countOfReadsProcessed + " reads.");

            int indexOfFirstComma = line.indexOf(',');
            String idOfRead = line.substring(1);
            if (indexOfFirstComma > 0)
                idOfRead = line.substring(1, indexOfFirstComma);


            ExtendedReadMappingSetForRead setOfSortedExtendedReadMappings = new ExtendedReadMappingSetForRead(idOfRead, new ExtendedReadMappingByScoreComparator());

            ExtendedReadMapping readMappings[] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(line.substring(1));

            String sequenceOfRead = bufferedRandomAccessFileMAX.readNextLine();
            line = bufferedRandomAccessFileMAX.readNextLine();


            int countOfAlignmentsForThisReadWithMinAlignScore = 0;

            for (int indexOfMapping = 0; indexOfMapping < readMappings.length; indexOfMapping++) {

                if (readMappings[indexOfMapping].getScore() >= minAlignmentScoreForReportingAlignment) {
                    countOfAlignmentsForThisReadWithMinAlignScore++;
                    setOfSortedExtendedReadMappings.addMapping(readMappings[indexOfMapping]);
                }
            }

            if (countOfAlignmentsForThisReadWithMinAlignScore >= minAlignmentsWithMinScoreRequiredBeforeReportingRead
                    && countOfAlignmentsForThisReadWithMinAlignScore <= maxAlignmentsWithMinScoreAllowedBeforeNotReportingRead) {

                countOfReadsWithNumberOfAlignmentsInRequiredRange++;

                Iterator<ExtendedReadMapping> iteratorOverExtendedReadMappings = setOfSortedExtendedReadMappings.getSetOfSortedReadMappings().iterator();
                ExtendedReadMapping extendedReadMappingBestFound = null;
                ExtendedReadMapping extendedReadMappingSecondBestFound = null;

                if (iteratorOverExtendedReadMappings.hasNext())
                    extendedReadMappingBestFound = iteratorOverExtendedReadMappings.next();

                if (iteratorOverExtendedReadMappings.hasNext())
                    extendedReadMappingSecondBestFound = iteratorOverExtendedReadMappings.next();


                if (extendedReadMappingSecondBestFound == null || extendedReadMappingBestFound.getScore() - extendedReadMappingSecondBestFound.getScore() >= minScoreGapToSecondBestAlignment) {

                    //ReadMapping readMapping = new ReadMapping(extendedReadMappingBestFound.getIndexOfMatchingReferenceSequence(), extendedReadMappingBestFound.getPositionOfAlignmentStartInReferenceSequence());
                    //readMapping.setIdOfRead(extendedReadMappingBestFound.getIdOfRead());
                    Read read = new Read(extendedReadMappingBestFound.getIdOfRead() + extendedReadMappingBestFound.toString(), sequenceOfRead);
                    mapSortedReadMappingToMergedMaxFileEntry.put(extendedReadMappingBestFound, read);

                }

            }

        }

        bufferedRandomAccessFileMAX.close();

        BufferedWriter writerMAXSortedOutput = new BufferedWriter(new FileWriter(fileMAXSortedOutput));
        Iterator<ExtendedReadMapping> iteratorOverSortedReadMappings = mapSortedReadMappingToMergedMaxFileEntry.keySet().iterator();
        while (iteratorOverSortedReadMappings.hasNext()) {
            ExtendedReadMapping readMapping = iteratorOverSortedReadMappings.next();
            Read read = mapSortedReadMappingToMergedMaxFileEntry.get(readMapping);
            writerMAXSortedOutput.write(">" + read.getNameOfRead());
            writerMAXSortedOutput.newLine();
            writerMAXSortedOutput.write(read.getSequenceOfRead());
            writerMAXSortedOutput.newLine();
        }
        writerMAXSortedOutput.close();

        System.out.println("countOfReadsProcessed:\t" + countOfReadsProcessed);
        System.out.println("countOfReadsWithNumberOfAlignmentsInRequiredRange:\t" + countOfReadsWithNumberOfAlignmentsInRequiredRange);


    }

}
