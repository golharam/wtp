package com.lifetechnologies.solid.wt;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeatureList;
import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.ContiguousGenomicFeatureComparator;
import com.lifetechnologies.solid.wt.Strand;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.TreeSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * User: tuchbb
 * Date: Mar 3, 2009
 * Time: 9:11:23 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExonList extends ContiguousGenomicFeatureList {

    public ExonList(File fileRefGeneTable, int sizeOfBins) throws Exception {
        super(sizeOfBins);

        // load all the features into a temporary data structure
        HashMap<String, TreeSet<ContiguousGenomicFeature>> mapContigNameToSortedExonSet = new HashMap<String, TreeSet<ContiguousGenomicFeature>>();
        HashMap<String, Integer> mapContigNameToMaxCoordinate = new HashMap<String, Integer>();

        BufferedReader readerRefGeneTable = new BufferedReader(new FileReader(fileRefGeneTable));
        String line = null;
        while ((line = readerRefGeneTable.readLine()) != null) {
            if (!line.startsWith("#") && line.length() > 0) {
                String tokens[] = line.split("\t");

                String idRefSeq = tokens[1];
                String nameOfGene = tokens[12];
                String nameOfContig = tokens[2];
                Strand strand = Strand.toStrand(tokens[3].toCharArray()[0]);

                String stringAdditionalInfo = "ID=" + nameOfGene + ";refSeq=" + idRefSeq;

                int coordinateCDSStart = Integer.parseInt(tokens[6]) +1;    // convert 0-based to 1-based
                int coordinateCDSEnd = Integer.parseInt(tokens[7]);         // refGene ends are non-inclusive

                String tokensExonStarts[] = tokens[9].split(",");
                String tokensExonEnds[] = tokens[10].split(",");
                String tokensExonFrames[] = tokens[15].split(",");

                if (tokensExonStarts.length != tokensExonFrames.length) {
                    System.out.println("Skipping " + idRefSeq + " because # of exon starts (" + tokensExonStarts.length + ") does not equal # of frames listed (" + tokensExonFrames.length + ")");
                    continue;
                }

                Exon exonOneBefore = null;
                Exon exonCurrent = null;

                for (int i = 0; i < tokensExonStarts.length; i++) {
                    int coordinateStartExon = Integer.parseInt(tokensExonStarts[i]) +1; // convert 0-based to 1-based
                    int coordinateEndExon = Integer.parseInt(tokensExonEnds[i]);    // refGene ends are non-inclusive

                    exonCurrent = new Exon(nameOfContig, strand, coordinateStartExon, coordinateEndExon);

                    Frame frameOfCodingSequenceForFeature = Frame.toFrame(Integer.parseInt(tokensExonFrames[i]));
                    exonCurrent.setFrameOfCodingSequence(frameOfCodingSequenceForFeature);
                    if (frameOfCodingSequenceForFeature != Frame.NONE) {

                        int positionCodingStartInExon = 1;
                        if (coordinateCDSStart > coordinateStartExon)
                            positionCodingStartInExon += coordinateCDSStart - coordinateStartExon;
                        exonCurrent.setPositionOfCodingSequenceStart(positionCodingStartInExon);

                        int positionCodingEndInExon = exonCurrent.getLength();
                        if (coordinateCDSEnd < coordinateEndExon)
                            positionCodingEndInExon -= coordinateEndExon - coordinateCDSEnd;
                        exonCurrent.setPositionOfCodingSequenceEnd(positionCodingEndInExon);
                    }

                    if (strand == Strand.POSITIVE) {
                        exonCurrent.setPreviousExon(exonOneBefore);
                        if (exonOneBefore!= null)    exonOneBefore.setNextExon(exonCurrent);
                    } else {
                        exonCurrent.setNextExon(exonOneBefore);
                        if (exonOneBefore!= null)    exonOneBefore.setPreviousExon(exonCurrent);
                    }
                    exonCurrent.setAdditionalInfo(stringAdditionalInfo);

                    if (mapContigNameToSortedExonSet.containsKey(nameOfContig))
                        mapContigNameToSortedExonSet.get(nameOfContig).add(exonCurrent);
                    else {
                        TreeSet<ContiguousGenomicFeature> setOfFeatures = new TreeSet<ContiguousGenomicFeature>(new ContiguousGenomicFeatureComparator());
                        setOfFeatures.add(exonCurrent);
                        mapContigNameToSortedExonSet.put(nameOfContig, setOfFeatures);
                    }
                    if (mapContigNameToMaxCoordinate.containsKey(nameOfContig)) {
                        if (mapContigNameToMaxCoordinate.get(nameOfContig) < coordinateEndExon)
                            mapContigNameToMaxCoordinate.put(nameOfContig, coordinateEndExon);
                    } else
                        mapContigNameToMaxCoordinate.put(nameOfContig, coordinateEndExon);

                    exonOneBefore = exonCurrent;
                }

            }
        }
        readerRefGeneTable.close();

        super.loadDataStructure(sizeOfBins, mapContigNameToSortedExonSet, mapContigNameToMaxCoordinate);

    }

    public static void main(String args[]) throws Exception  {

        ContiguousGenomicFeatureList listOfExons = new ExonList(new File("/home/tuchbb/data/species/h_sapiens/ucsc_hg18_081115/refGene.081115.filtered.first100Entries"), 1000);
        HashSet<ContiguousGenomicFeature> featuresOverlapping = listOfExons.getSetOfFeaturesOverlapping("chr1", 869500, 869600).getFeaturesThatOverlapThisOneWithOverlapAnnotated();
        //ContiguousGenomicFeatureList listOfGenes = new ContiguousGenomicFeatureList(new File("/home/tuchbb/data/species/h_sapiens/ucsc_hg18_081115/gene_annotation.filtered.081115.gff"), 1000, 0, true, true, false, true);
        //HashSet<ContiguousGenomicFeature> featuresOverlapping = listOfGenes.getSetOfFeaturesOverlapping("chr5", 132000000, 133000000).getFeaturesThatOverlapThisOneWithOverlapAnnotated();
        Iterator<ContiguousGenomicFeature> iterator = featuresOverlapping.iterator();
        while (iterator.hasNext()) {
            System.out.println(iterator.next());
        }
    }


}
