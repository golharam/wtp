package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Iterator;
import java.util.HashMap;
import java.util.TreeMap;

/**
 * User: tuchbb
 * Date: Mar 3, 2009
 * Time: 11:58:03 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class GenomePositionAnnotator {


    public static void main(String args[]) throws Exception  {

        BufferedReader readerPositionsToAnnotate = new BufferedReader(new FileReader(args[0]));
        ExonList listOfExons = new ExonList(new File(args[1]), 1000);
        CodonTable codonTable = new CodonTable(new File(args[2]));
        FastaDatabase fastaDatabaseReference = new FastaDatabase(new File(args[3]));

        char charNucleotides[] = { 'A', 'C', 'G', 'T'};

        TreeMap<ContiguousGenomicFeature, String> mapFeaturePositionToLine = new TreeMap<ContiguousGenomicFeature, String>(new ContiguousGenomicFeatureComparator());
        String line;
        while ((line = readerPositionsToAnnotate.readLine()) != null) {
            String tokens[] = line.split("\t");
            String nameOfContig = tokens[0];
            int coordinate = Integer.parseInt(tokens[1].replaceAll(",", ""));
            mapFeaturePositionToLine.put(new ContiguousGenomicFeature(nameOfContig, Strand.EITHER, coordinate, coordinate), line);
        }
        readerPositionsToAnnotate.close();

        StringBuffer sequenceOfContig = null;
        String nameOfLastContig = "";
        Iterator<ContiguousGenomicFeature> iteratorOverCoordiantes = mapFeaturePositionToLine.keySet().iterator();
        while (iteratorOverCoordiantes.hasNext()) {
            ContiguousGenomicFeature featureCoordinate = iteratorOverCoordiantes.next();
            line = mapFeaturePositionToLine.get(featureCoordinate);
            //while ((line = readerPositionsToAnnotate.readLine()) != null) {
            String tokens[] = line.split("\t");
            String nameOfContig = tokens[0];
            
            if (!nameOfContig.equals(nameOfLastContig))     // FIXME: If input file is not sorted by chromosomal coordinate then this will be very slow!
                sequenceOfContig = fastaDatabaseReference.getSequenceByHeader(nameOfContig);
            nameOfLastContig = nameOfContig;

            int coordinate = Integer.parseInt(tokens[1].replaceAll(",", ""));
            HashSet<String> setOfAnnotationsForPosition = new HashSet<String>();

            HashSet<ContiguousGenomicFeature> setOfOverlappingExons = listOfExons.getSetOfFeaturesOverlapping(nameOfContig, coordinate, coordinate).getFeaturesThatOverlapThisOne();
            Iterator<ContiguousGenomicFeature> iteratorOverOverlappingExons = setOfOverlappingExons.iterator();
            while (iteratorOverOverlappingExons.hasNext()) {
                Exon exonOverlapping = (Exon)iteratorOverOverlappingExons.next();
                String annotation = null;
                if (exonOverlapping.getFrameOfCodingSequence() == Frame.NONE) {

                    boolean hitCodingExonFivePrime = false;
                    Exon previousExon = exonOverlapping;
                    while ((previousExon = previousExon.getPreviousExon()) != null && !hitCodingExonFivePrime) {
                        if (previousExon.getFrameOfCodingSequence() != Frame.NONE)
                            hitCodingExonFivePrime = true;
                    }
                    boolean hitCodingExonThreePrime = false;
                    Exon nextExon = exonOverlapping;
                    while ((nextExon = nextExon.getNextExon()) != null && !hitCodingExonThreePrime) {
                        if (nextExon.getFrameOfCodingSequence() != Frame.NONE)
                            hitCodingExonThreePrime = true;
                    }

                    if (!hitCodingExonFivePrime && !hitCodingExonThreePrime)
                        annotation = "non-coding";
                    else if (!hitCodingExonFivePrime && hitCodingExonThreePrime)
                        annotation = "5'utr";
                    else if (hitCodingExonFivePrime && !hitCodingExonThreePrime)
                        annotation = "3'utr";

                } else {

                    if (exonOverlapping.getStrand() == Strand.POSITIVE) {
                        if (coordinate < (exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceStart() -1))
                            annotation = "5'utr";
                        else if (coordinate > (exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceEnd() -1))
                            annotation = "3'utr";
                        else {  // extract the relevant codon and determine if synonymous or non-synonymous

                            int coordinateOfCodingSequenceStart = exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceStart() -1;
                            Frame frameOfCoordinateInExon = Frame.toFrame((coordinate - coordinateOfCodingSequenceStart + exonOverlapping.getFrameOfCodingSequence().toInt()) % 3);
                            String codon = sequenceOfContig.substring(coordinate -1, coordinate).toUpperCase();
                            HashMap<Character, String> mapNucleotideChangeToCodon = new HashMap<Character, String>();
                            for (int i = 0; i < charNucleotides.length; i++)
                                mapNucleotideChangeToCodon.put(charNucleotides[i], charNucleotides[i] + "");

                            int coordinateLocal = coordinate;
                            Exon exonLocal = exonOverlapping;
                            for (int indexLowerFrames = 0; indexLowerFrames < frameOfCoordinateInExon.toInt(); indexLowerFrames++) {
                                coordinateLocal--;
                                if (coordinateLocal < exonLocal.getCoordinateOfStart()) {
                                    exonLocal = exonLocal.getPreviousExon();
                                    coordinateLocal = exonLocal.getCoordinateOfEnd();
                                }
                                String nucleotideToAdd = sequenceOfContig.substring(coordinateLocal -1, coordinateLocal).toLowerCase();
                                codon = nucleotideToAdd + codon;
                                for (int i = 0; i < charNucleotides.length; i++)
                                    mapNucleotideChangeToCodon.put(charNucleotides[i], nucleotideToAdd + mapNucleotideChangeToCodon.get(charNucleotides[i]));
                            }
                            exonLocal = exonOverlapping;
                            coordinateLocal = coordinate;
                            for (int indexHigherFrames = 0; indexHigherFrames < 2 - frameOfCoordinateInExon.toInt(); indexHigherFrames++) {
                                coordinateLocal++;
                                if (coordinateLocal > exonLocal.getCoordinateOfEnd()) {
                                    exonLocal = exonLocal.getNextExon();
                                    coordinateLocal = exonLocal.getCoordinateOfStart();
                                }
                                String nucleotideToAdd = sequenceOfContig.substring(coordinateLocal -1, coordinateLocal).toLowerCase();
                                codon += nucleotideToAdd;
                                for (int i = 0; i < charNucleotides.length; i++)
                                    mapNucleotideChangeToCodon.put(charNucleotides[i], mapNucleotideChangeToCodon.get(charNucleotides[i]) + nucleotideToAdd);
                            }

                            String aminoAcidReference = codonTable.translateCodon(codon);
                            annotation = "coding" + "_" + codon;
                            boolean foundPossibleAAChange = false;
                            for (int i = 0; i < charNucleotides.length; i++) {
                                String aminoAcidCurrent = codonTable.translateCodon(mapNucleotideChangeToCodon.get(charNucleotides[i]));
                                if (!aminoAcidReference.equalsIgnoreCase(aminoAcidCurrent)) {
                                    foundPossibleAAChange = true;
                                    annotation += "," + mapNucleotideChangeToCodon.get(charNucleotides[i]) + ":" + aminoAcidReference + "->" + aminoAcidCurrent;
                                }
                            }
                            if (!foundPossibleAAChange)
                                annotation += "-synonymous";

                        }

                    } else if (exonOverlapping.getStrand() == Strand.NEGATIVE) {

                        if (coordinate < (exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceStart() -1))
                            annotation = "3'utr";
                        else if (coordinate > (exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceEnd() -1))
                            annotation = "5'utr";
                        else {  // extract the relevant codon and determine if synonymous or non-synonymous

                            int coordinateOfCodingSequenceEnd = exonOverlapping.getCoordinateOfStart() + exonOverlapping.getPositionOfCodingSequenceEnd() -1;
                            Frame frameOfCoordinateInExon = Frame.toFrame((coordinateOfCodingSequenceEnd - coordinate + exonOverlapping.getFrameOfCodingSequence().toInt()) % 3);
                            String codon = SequenceUtilities.reverseComplement(sequenceOfContig.substring(coordinate -1, coordinate)).toUpperCase();
                            HashMap<Character, String> mapNucleotideChangeToCodon = new HashMap<Character, String>();
                            for (int i = 0; i < charNucleotides.length; i++)
                                mapNucleotideChangeToCodon.put(charNucleotides[i], charNucleotides[i] + "");

                            int coordinateLocal = coordinate;
                            Exon exonLocal = exonOverlapping;
                            for (int indexLowerFrames = 0; indexLowerFrames < frameOfCoordinateInExon.toInt(); indexLowerFrames++) {
                                coordinateLocal++;
                                if (coordinateLocal > exonLocal.getCoordinateOfEnd()) {
                                    exonLocal = exonLocal.getPreviousExon();
                                    coordinateLocal = exonLocal.getCoordinateOfStart();
                                }
                                String nucleotideToAdd = SequenceUtilities.reverseComplement(sequenceOfContig.substring(coordinateLocal -1, coordinateLocal));
                                codon = nucleotideToAdd + codon;
                                for (int i = 0; i < charNucleotides.length; i++)
                                    mapNucleotideChangeToCodon.put(charNucleotides[i], nucleotideToAdd + mapNucleotideChangeToCodon.get(charNucleotides[i]));
                            }
                            exonLocal = exonOverlapping;
                            coordinateLocal = coordinate;
                            for (int indexHigherFrames = 0; indexHigherFrames < 2 - frameOfCoordinateInExon.toInt(); indexHigherFrames++) {
                                coordinateLocal--;
                                if (coordinateLocal < exonLocal.getCoordinateOfStart()) {
                                    exonLocal = exonLocal.getNextExon();
                                    coordinateLocal = exonLocal.getCoordinateOfEnd();
                                }
                                String nucleotideToAdd = SequenceUtilities.reverseComplement(sequenceOfContig.substring(coordinateLocal -1, coordinateLocal));
                                codon += nucleotideToAdd;
                                for (int i = 0; i < charNucleotides.length; i++)
                                    mapNucleotideChangeToCodon.put(charNucleotides[i], mapNucleotideChangeToCodon.get(charNucleotides[i]) + nucleotideToAdd);
                            }

                            String aminoAcidReference = codonTable.translateCodon(codon);
                            annotation = "coding" + "_" + codon;
                            boolean foundPossibleAAChange = false;
                            for (int i = 0; i < charNucleotides.length; i++) {
                                String aminoAcidCurrent = codonTable.translateCodon(mapNucleotideChangeToCodon.get(charNucleotides[i]));
                                if (!aminoAcidReference.equalsIgnoreCase(aminoAcidCurrent)) {
                                    foundPossibleAAChange = true;
                                    annotation += "," + mapNucleotideChangeToCodon.get(charNucleotides[i]) + ":" + aminoAcidReference + "->" + aminoAcidCurrent;
                                }
                            }
                            if (!foundPossibleAAChange)
                                annotation += "-synonymous";

                        }

                    }
                }

                if (annotation != null) {
                    String nameOfGene = exonOverlapping.getAdditionalInfo().split(";")[0].substring(3);
                    setOfAnnotationsForPosition.add(nameOfGene + "_" + annotation);
                }
            }

            // if the coordinate lies outside of any of the exons, determine if it is intronic or intergenic
            if (setOfOverlappingExons.size() == 0) {
                String annotation = null;
                HashSet<ContiguousGenomicFeature> setOfFlankingExons = listOfExons.getSetOfFeaturesImmediatelyFlanking(nameOfContig, coordinate, coordinate);
                Iterator<ContiguousGenomicFeature> iteratorOverFlankingExons = setOfFlankingExons.iterator();
                while (iteratorOverFlankingExons.hasNext()) {
                    Exon exonFlanking = (Exon)iteratorOverFlankingExons.next();
                    if (coordinate > exonFlanking.getCoordinateOfEnd()) {

                        if (exonFlanking.getStrand() == Strand.POSITIVE && exonFlanking.getNextExon() == null)
                            annotation = "intergenic";
                        else if (exonFlanking.getStrand() == Strand.NEGATIVE && exonFlanking.getPreviousExon() == null)
                            annotation = "intergenic";
                        else {
                            annotation = "intronic";
                            if (Math.abs(coordinate - exonFlanking.getCoordinateOfEnd()) < 10)
                                annotation += "-nearSpliceSite";
                        }
                    } else if (coordinate < exonFlanking.getCoordinateOfStart()) {

                        if (exonFlanking.getStrand() == Strand.POSITIVE && exonFlanking.getPreviousExon() == null)
                            annotation = "intergenic";
                        else if (exonFlanking.getStrand() == Strand.NEGATIVE && exonFlanking.getNextExon() == null)
                            annotation = "intergenic";
                        else {
                            annotation = "intronic";
                            if (Math.abs(coordinate - exonFlanking.getCoordinateOfStart()) < 10)
                                annotation += "-nearSpliceSite";
                        }
                    }

                    if (annotation != null) {
                        String nameOfGene = exonFlanking.getAdditionalInfo().split(";")[0].substring(3);
                        setOfAnnotationsForPosition.add(nameOfGene + "_" + annotation);
                    }
                }
            }

            System.out.print(nameOfContig + "\t" + coordinate + "\t");
            Iterator<String> iteratorOverAnnotations = setOfAnnotationsForPosition.iterator();
            while (iteratorOverAnnotations.hasNext())
                System.out.print(iteratorOverAnnotations.next() + ";");
            System.out.println();

        }
        //readerPositionsToAnnotate.close();


    }
}
