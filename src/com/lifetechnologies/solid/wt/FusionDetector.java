package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.HashMap;

/**
 * User: tuchbb
 * Date: Aug 21, 2008
 * Time: 9:36:16 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class FusionDetector {


    static void detectFusions(File fileInputSummaryA, File fileInputSummaryB, File fileOutputFusionDetection, boolean verbose) throws IOException {

        File fileInputSummarySmaller = fileInputSummaryA;
        File fileInputSummaryLarger = fileInputSummaryB;
        boolean flippedFiles = false;
        if (fileInputSummarySmaller.length() > fileInputSummaryLarger.length()) {
            fileInputSummarySmaller = fileInputSummaryB;
            fileInputSummaryLarger = fileInputSummaryA;
            flippedFiles = true;
        }

        HashMap<String, ReadMapping> mapReadIdToReadMapping = new HashMap<String, ReadMapping>();
        BufferedReader readerOutputSummarySmaller = new BufferedReader(new FileReader(fileInputSummarySmaller));
        String line = readerOutputSummarySmaller.readLine();  // drop header
        while ((line = readerOutputSummarySmaller.readLine()) != null) {
            String tokens[] = line.split("\t");
            // FIXME: if a read maps to multiple targets, then this will currently keep only one such mapping
            if (mapReadIdToReadMapping.containsKey(tokens[0]))
                System.err.println("warning: multiple matches for same read id; previous match will be replaced");
            mapReadIdToReadMapping.put(tokens[0], new ReadMapping(tokens[1], Integer.parseInt(tokens[2]), Integer.parseInt(tokens[3])));
        }
        readerOutputSummarySmaller.close();

        BufferedWriter writerOutputFusionDetection = new BufferedWriter(new FileWriter(fileOutputFusionDetection));
        writerOutputFusionDetection.write("idOfMatchingReferenceSequenceA\tpositionOfMatchInReferenceSequenceA\tcountMismatchesA\tidOfMatchingReferenceSequenceB\tpositionOfMatchInReferenceSequenceB\tcountMismatchesB");
        writerOutputFusionDetection.newLine();
        int countReadMatchesSameGene = 0;
        int countReadMatchesSameGeneButDifferentPosition = 0;
        int countReadMatchesDifferentGene = 0;
        BufferedReader readerOutputSummaryLarger = new BufferedReader(new FileReader(fileInputSummaryLarger));
        line = readerOutputSummaryLarger.readLine(); // drop header
        while ((line = readerOutputSummaryLarger.readLine()) != null) {
            String tokens[] = line.split("\t");
            if (mapReadIdToReadMapping.containsKey(tokens[0])) {
                ReadMapping readMapping = mapReadIdToReadMapping.get(tokens[0]);
                if (readMapping.getIdOfMatchingReferenceSequence().equalsIgnoreCase(tokens[1])) {
                    countReadMatchesSameGene++;
                    if (readMapping.getPositionOfAlignmentStartInReferenceSequence() != Integer.parseInt(tokens[2]))
                        countReadMatchesSameGeneButDifferentPosition++;
                } else {
                    countReadMatchesDifferentGene++;
                    if (flippedFiles)   writerOutputFusionDetection.write(tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + readMapping.getIdOfMatchingReferenceSequence() + "\t" + readMapping.getPositionOfAlignmentStartInReferenceSequence() + "\t" + readMapping.getNumberOfMismatches());
                    else    writerOutputFusionDetection.write(readMapping.getIdOfMatchingReferenceSequence() + "\t" + readMapping.getPositionOfAlignmentStartInReferenceSequence() + "\t" + readMapping.getNumberOfMismatches() + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3]);
                    writerOutputFusionDetection.newLine();
                }
            }
        }
        readerOutputSummaryLarger.close();
        writerOutputFusionDetection.close();

        if (verbose) {
            System.out.println("countReadMatchesSameGene:\t" + countReadMatchesSameGene);
            System.out.println("countReadMatchesSameGeneButDifferentPosition:\t" + countReadMatchesSameGeneButDifferentPosition);
            System.out.println("countReadMatchesDifferentGene:\t" + countReadMatchesDifferentGene);
        }
    }

    @SuppressWarnings("unused")
    public static void main(String args[]) throws IOException {

        //File fileTrascriptIds = new File("/home/tuchbb/tuchbb_on_vast/projects/GEx/translocations/fusion_transcript_ids.txt");
        File fileFusionPoints = new File("/home/tuchbb/tuchbb_on_vast/projects/GEx/translocations/translocation_190assays_info.txt");
        //File fileReferenceFasta = new File("/home/tuchbb/data/species/h_sapiens/Homo_sapiens.NCBI36.38.apr.cdna.fa");
        File fileReferenceFasta = new File("/home/tuchbb/tuchbb_on_vast/projects/GEx/translocations/fusion_transcripts.from190Assays.ffn");
        //File fileReferenceFastaFiltered = new File("/home/tuchbb/data/species/h_sapiens/Homo_sapiens.NCBI36.38.apr.cdna.justKnownFusionTranscripts.fa");
        File fileReferenceFastaFiltered = new File("/home/tuchbb/tuchbb_on_vast/projects/GEx/translocations/fusion_transcripts.from190Assays.justPaddedFusionRegion.ffn");

        FastaDatabase fastaDBReference = new FastaDatabase(fileReferenceFasta);

        //FIXME: move to com.lifetechnologies.solid.wt.FastaDatabase class
        // com.lifetechnologies.solid.wt.FastaDatabase fastaDBReference = extractSubsetOfFastaDBtoNewFastaDBFile(fileTrascriptIds, fastaDBReference, fileReferenceFastaFiltered);

        //RichSequence sequence = fastaDBReference.getSequenceByHeaderPrefix("ENST00000305877");
        //SymbolList symL = RNATools.createRNA("auggcaccguccagauu");
        //System.out.println(sequence.subStr(2,50));

        
    }
     
}
