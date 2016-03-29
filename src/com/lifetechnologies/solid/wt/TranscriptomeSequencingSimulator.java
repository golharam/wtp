package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;

/**
 * User: tuchbb
 * Date: Sep 22, 2008
 * Time: 1:32:43 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 *
 * FIXME: UNTESTED!!!!
 */
public class TranscriptomeSequencingSimulator {

    private static File FILE_TRANSCRIPTS_FASTA = new File("/home/tuchbb/data/species/h_sapiens/ucsc_hg18/refSeq_mRNA.fa");
//    private static File FILE_POSITION_SPECIFIC_ERROR_PROFILE = new File("/home/tuchbb/projects/simulatedReads/errorProfile.cloonan.txt");
    private static File FILE_SIMULATED_READS_CSFASTA = new File("/home/tuchbb/data/species/h_sapiens/simulatedReads/simulatedReads.ucscHG18_refSeqMRNA.fa");

    private static int DYNAMIC_RANGE_IN_LOG_TEN = 8;
    private static int LENGTH_OF_READ = 50;
    private static int NUMBER_OF_READS_TO_GENERATE = 200000000;

    private static char NUCLEOTIDE_AT_FIVE_PRIME_END_OF_READ = 'T';

    public static void main(String args[]) throws Exception {

        FastaDatabase fastaDatabaseOfTranscipts = new FastaDatabase(FILE_TRANSCRIPTS_FASTA);
        HashMap<String, String> mapHeaderToSequence = fastaDatabaseOfTranscipts.getHeaderToSequenceMap(true);
        String headers[] = (String[])mapHeaderToSequence.keySet().toArray();
        double maxExpressionLevel = Math.pow(10, DYNAMIC_RANGE_IN_LOG_TEN);
        double probabilitiesOfTranscriptSampling[] = new double[headers.length];
        for (int i = 0; i < probabilitiesOfTranscriptSampling.length; i++)
            probabilitiesOfTranscriptSampling[i] = Math.pow(10, Math.random() * DYNAMIC_RANGE_IN_LOG_TEN) / maxExpressionLevel;

        BufferedWriter writerSimulatedReadsFile = new BufferedWriter(new FileWriter(FILE_SIMULATED_READS_CSFASTA));
        int countNumberOfReadsGenerated = 0;
        while (countNumberOfReadsGenerated < NUMBER_OF_READS_TO_GENERATE) {
            int indexOfRandomTranscript = (int) (Math.random() * probabilitiesOfTranscriptSampling.length);
            if (Math.random() > probabilitiesOfTranscriptSampling[indexOfRandomTranscript]) {
                String sequenceOfTranscriptInBases = mapHeaderToSequence.get(headers[indexOfRandomTranscript]);
                int indexOfRandomPositionInTranscript = (int) (Math.random() * (sequenceOfTranscriptInBases.length() - LENGTH_OF_READ +1));
                String sequenceOfReadInBases = NUCLEOTIDE_AT_FIVE_PRIME_END_OF_READ + sequenceOfTranscriptInBases.substring(indexOfRandomPositionInTranscript, indexOfRandomPositionInTranscript + LENGTH_OF_READ);
                String sequenceOfReadInColors = BaseSpaceToColorSpaceTranslator.translateSequence(sequenceOfReadInBases);

                char sequenceOfReadInColorsAsArray[] = sequenceOfReadInColors.toCharArray();
                for (int i = 0; i < sequenceOfReadInColorsAsArray.length; i++) {
                    // FIXME: add errors and SNPs

                }

                writerSimulatedReadsFile.write(">0_0_" + headers[indexOfRandomTranscript] + "_" + indexOfRandomPositionInTranscript);
                writerSimulatedReadsFile.newLine();
                writerSimulatedReadsFile.write(NUCLEOTIDE_AT_FIVE_PRIME_END_OF_READ + new String(sequenceOfReadInColorsAsArray));
                writerSimulatedReadsFile.newLine();

                countNumberOfReadsGenerated++;
            }
        }
        writerSimulatedReadsFile.close();

    }
}
