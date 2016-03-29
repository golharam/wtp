package com.lifetechnologies.solid.wt;


import java.util.*;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;

/**
 * User: tuchbb
 * Date: Jan 15, 2009
 * Time: 11:14:15 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ContiguousGenomicFeatureList {

    // the genome is binned with size = 1000, and features with any coordinates overlapping that bin are stored in this structure
    private HashMap<String, HashSet<ContiguousGenomicFeature>[]> mapContigNameToArrayOfBinnedFeatureSets;
    private TreeSet<ContiguousGenomicFeature> setOfAllFeatures;
    private int sizeOfBins;

    public ContiguousGenomicFeatureList(int sizeOfBins) {
        this.sizeOfBins = sizeOfBins;
        this.setOfAllFeatures = new TreeSet<ContiguousGenomicFeature>(new ContiguousGenomicFeatureComparator());        
    }

    public ContiguousGenomicFeatureList(File fileAnnotationGFF, int sizeOfBins, int numberOfBasesToLengthenFeaturesBy,
                                        boolean storeSourceField, boolean storeTypeField, boolean storeScoreField, boolean storeAttributesField) throws Exception {

        this.sizeOfBins = sizeOfBins;
        this.setOfAllFeatures = new TreeSet<ContiguousGenomicFeature>(new ContiguousGenomicFeatureComparator());

        // load all the features into a temporary data structure
        HashMap<String, TreeSet<ContiguousGenomicFeature>> mapContigNameToSortedFeatueSet = new HashMap<String, TreeSet<ContiguousGenomicFeature>>();
        HashMap<String, Integer> mapContigNameToMaxCoordinate = new HashMap<String, Integer>();

        BufferedReader readerAnnotationGFF = new BufferedReader(new FileReader(fileAnnotationGFF));
        String line = null;
        while ((line = readerAnnotationGFF.readLine()) != null) {
            if (!line.startsWith("#") && line.length() > 0) {
                String tokens[] = line.split("\t");
                Strand strand = Strand.toStrand(tokens[6].toCharArray()[0]);

                //if (strand != null) {
                String nameOfContig = tokens[0];
                int coordinateStart = Integer.parseInt(tokens[3]);
                coordinateStart = Math.max(1, coordinateStart - numberOfBasesToLengthenFeaturesBy);
                int coordinateEnd = Integer.parseInt(tokens[4]) + numberOfBasesToLengthenFeaturesBy;
                ContiguousGenomicFeature feature = new ContiguousGenomicFeature(nameOfContig, strand, coordinateStart, coordinateEnd);
                if (storeSourceField && !tokens[1].equals("."))   feature.setLabel(tokens[1]);
                if (storeTypeField && !tokens[2].equals(".")) feature.setType(tokens[2]);
                if (storeScoreField && !tokens[5].equals("."))    feature.setValueA(Double.parseDouble(tokens[5]));
                if (storeAttributesField && !tokens[8].equals("."))   feature.setAdditionalInfo(tokens[8]);
                if (mapContigNameToSortedFeatueSet.containsKey(nameOfContig))
                    mapContigNameToSortedFeatueSet.get(nameOfContig).add(feature);
                else {
                    TreeSet<ContiguousGenomicFeature> setOfFeatures = new TreeSet<ContiguousGenomicFeature>(new ContiguousGenomicFeatureComparator());
                    setOfFeatures.add(feature);
                    mapContigNameToSortedFeatueSet.put(nameOfContig, setOfFeatures);
                }
                if (mapContigNameToMaxCoordinate.containsKey(nameOfContig)) {
                    if (mapContigNameToMaxCoordinate.get(nameOfContig) < coordinateEnd)
                        mapContigNameToMaxCoordinate.put(nameOfContig, coordinateEnd);
                } else
                    mapContigNameToMaxCoordinate.put(nameOfContig, coordinateEnd);
                //}

            }
        }
        readerAnnotationGFF.close();

        loadDataStructure(sizeOfBins, mapContigNameToSortedFeatueSet, mapContigNameToMaxCoordinate);


    }


    protected void loadDataStructure(int sizeOfBins, HashMap<String, TreeSet<ContiguousGenomicFeature>> mapContigNameToSortedFeatueSet, HashMap<String, Integer> mapContigNameToMaxCoordinate) {
        // now that all the features are loaded, we can move them to a data structure which will be faster to search
        // the genome is binned with size = 1000, and features with any coordinates overlapping a given bin are referenced from that bin
        this.mapContigNameToArrayOfBinnedFeatureSets = new HashMap<String, HashSet<ContiguousGenomicFeature>[]>();
        Iterator<String> iteratorOverContigNames = mapContigNameToSortedFeatueSet.keySet().iterator();
        while (iteratorOverContigNames.hasNext()) {
            String nameOfContig = iteratorOverContigNames.next();
            int maxCoordianateWithFeatureInContig = mapContigNameToMaxCoordinate.get(nameOfContig);
            int numberOfBins = (int)Math.ceil((maxCoordianateWithFeatureInContig +1) / (double)sizeOfBins);
            @SuppressWarnings("unchecked")
            HashSet<ContiguousGenomicFeature> setsOfFeaturesByGenomicBin[] = (HashSet<ContiguousGenomicFeature>[])new HashSet[numberOfBins];
            TreeSet<ContiguousGenomicFeature> sortedSetOfFeatures = mapContigNameToSortedFeatueSet.get(nameOfContig);
            Iterator<ContiguousGenomicFeature> iteratorOverFeatures = sortedSetOfFeatures.iterator();
            while (iteratorOverFeatures.hasNext()) {
                ContiguousGenomicFeature feature = iteratorOverFeatures.next();
                //System.out.println(nameOfContig + "\t" + maxCoordianateWithFeatureInContig + "\t" + numberOfBins + "\t" + feature);
                this.setOfAllFeatures.add(feature);
                int binFirst = feature.getCoordinateOfStart() / sizeOfBins;
                int binLast = feature.getCoordinateOfEnd() / sizeOfBins;
                for (int i = binFirst; i <= binLast ; i++) {
                    if (setsOfFeaturesByGenomicBin[i] == null)
                        setsOfFeaturesByGenomicBin[i] = new HashSet<ContiguousGenomicFeature>();
                    setsOfFeaturesByGenomicBin[i].add(feature);
                }
            }
            this.mapContigNameToArrayOfBinnedFeatureSets.put(nameOfContig, setsOfFeaturesByGenomicBin);
        }
    }


    public OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(ContiguousGenomicFeature feature, boolean requireSameStrand, double minOverlapFraction) throws CloneNotSupportedException {
        Strand strand = Strand.EITHER;
        if (requireSameStrand && feature.getStrand() != null)
            strand = feature.getStrand();
        return getSetOfFeaturesOverlapping(feature.getIdOfReferenceSequence(), strand, feature.getCoordinateOfStart(), feature.getCoordinateOfEnd(), minOverlapFraction, false);
    }

    public OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(ContiguousGenomicFeature feature, boolean requireSameStrand, boolean requireFeaturesToFullySpanThisRange) throws CloneNotSupportedException {
        Strand strand = Strand.EITHER;
        if (requireSameStrand && feature.getStrand() != null)
            strand = feature.getStrand();
        return getSetOfFeaturesOverlapping(feature.getIdOfReferenceSequence(), strand, feature.getCoordinateOfStart(), feature.getCoordinateOfEnd(), null, requireFeaturesToFullySpanThisRange);
    }

    public OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(String nameOfContig, Strand strand, int coordinateStart, int coordinateEnd, boolean requireFeaturesToFullySpanThisRange) throws CloneNotSupportedException {
            return getSetOfFeaturesOverlapping(nameOfContig, strand, coordinateStart, coordinateEnd, null, requireFeaturesToFullySpanThisRange);
    }

    public OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(String nameOfContig, int coordinateStart, int coordinateEnd) throws CloneNotSupportedException {
        return getSetOfFeaturesOverlapping(nameOfContig, Strand.EITHER, coordinateStart, coordinateEnd);
    }

    public OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(String nameOfContig, Strand strand, int coordinateStart, int coordinateEnd) throws CloneNotSupportedException {
        return getSetOfFeaturesOverlapping(nameOfContig, strand, coordinateStart, coordinateEnd, null, false);
    }

    private OverlappingGenomicFeaturesSet getSetOfFeaturesOverlapping(String nameOfContig, Strand strand, int coordinateStart, int coordinateEnd,
                                                                      Double minOverlapFraction, boolean requireFeaturesToFullySpanThisRange) throws CloneNotSupportedException {

        OverlappingGenomicFeaturesSet setOfFeaturesOverlapping = new OverlappingGenomicFeaturesSet(coordinateStart, coordinateEnd);
        if (setOfAllFeatures.size() == 0) return setOfFeaturesOverlapping;
        HashSet<ContiguousGenomicFeature> setOfFeatureAlreadyConsidered = new HashSet<ContiguousGenomicFeature>();

        if (this.mapContigNameToArrayOfBinnedFeatureSets.containsKey(nameOfContig)) {
            HashSet<ContiguousGenomicFeature> setsOfFeaturesByGenomicBin[] = this.mapContigNameToArrayOfBinnedFeatureSets.get(nameOfContig);
            int binFirst = coordinateStart / this.sizeOfBins;
            int binLast = coordinateEnd / this.sizeOfBins;
            for (int binCurrent = binFirst; binCurrent <= binLast ; binCurrent++) {
                if (binCurrent < setsOfFeaturesByGenomicBin.length && setsOfFeaturesByGenomicBin[binCurrent] != null) {
                    Iterator<ContiguousGenomicFeature> iteratorOverFeaturesInBin = setsOfFeaturesByGenomicBin[binCurrent].iterator();
                    while (iteratorOverFeaturesInBin.hasNext()) {
                        ContiguousGenomicFeature feature = iteratorOverFeaturesInBin.next();
                        if (!setOfFeatureAlreadyConsidered.contains(feature) &&
                                (strand == Strand.EITHER || feature.getStrand() == Strand.EITHER || feature.getStrand() == strand)) {

                            setOfFeatureAlreadyConsidered.add(feature);

                            int lengthOfOverlapingSequence = 1 + Math.min(feature.getCoordinateOfEnd(), coordinateEnd) - Math.max(feature.getCoordinateOfStart(), coordinateStart);
                            double overlapFromRegionPerspective = lengthOfOverlapingSequence / (double) (coordinateEnd - coordinateStart +1);
                            double overlapFromAnnotatedFeaturePerspective = lengthOfOverlapingSequence / (double) (feature.getCoordinateOfEnd() - feature.getCoordinateOfStart() +1);

                            ContiguousGenomicFeature featureAnnotatedWithOverlap = (ContiguousGenomicFeature) feature.clone();
                            featureAnnotatedWithOverlap.setOverlapFractionFromMyPerspective(overlapFromAnnotatedFeaturePerspective);
                            featureAnnotatedWithOverlap.setOverlapFractionFromItsPerspective(overlapFromRegionPerspective);

                            if (minOverlapFraction == null) {
                                if (!requireFeaturesToFullySpanThisRange && !(feature.getCoordinateOfEnd() < coordinateStart || feature.getCoordinateOfStart() > coordinateEnd))
                                    setOfFeaturesOverlapping.addOverlappingFeature(featureAnnotatedWithOverlap, feature);
                                else if (requireFeaturesToFullySpanThisRange && (feature.getCoordinateOfStart() <= coordinateStart && feature.getCoordinateOfEnd() >= coordinateEnd))
                                    setOfFeaturesOverlapping.addOverlappingFeature(featureAnnotatedWithOverlap, feature);
                            } else if (overlapFromRegionPerspective >= minOverlapFraction || overlapFromAnnotatedFeaturePerspective >= minOverlapFraction)
                                    setOfFeaturesOverlapping.addOverlappingFeature(featureAnnotatedWithOverlap, feature);

                        }
                    }
                }
            }
        }
        return setOfFeaturesOverlapping;
    }

    public HashSet<ContiguousGenomicFeature> getSetOfFeaturesImmediatelyFlanking(ContiguousGenomicFeature feature) {
        return getSetOfFeaturesImmediatelyFlanking(feature.getIdOfReferenceSequence(), feature.getCoordinateOfStart(), feature.getCoordinateOfEnd());
    }

    public HashSet<ContiguousGenomicFeature> getSetOfFeaturesImmediatelyFlanking(String nameOfContig, int coordinateStart, int coordinateEnd) {

        HashSet<ContiguousGenomicFeature> setOfFeaturesFlanking = getFirstFeatureWithEndCoordinateLessThan(nameOfContig, coordinateStart);
        setOfFeaturesFlanking.addAll(getFirstFeatureWithStartCoordinateGreaterThan(nameOfContig, coordinateEnd));

        return setOfFeaturesFlanking;
    }

    public HashSet<ContiguousGenomicFeature> getFirstFeatureWithEndCoordinateLessThan(String nameOfContig, int coordinate) {

        HashSet<ContiguousGenomicFeature> setOfFeatures = new HashSet<ContiguousGenomicFeature>();
        HashSet<ContiguousGenomicFeature> setsOfFeaturesByGenomicBin[] = this.mapContigNameToArrayOfBinnedFeatureSets.get(nameOfContig);
        if (setsOfFeaturesByGenomicBin != null) {
            int bin = Math.min(coordinate / this.sizeOfBins, setsOfFeaturesByGenomicBin.length -1);
            while (setOfFeatures.size() == 0 && bin >= 0) {
                if (setsOfFeaturesByGenomicBin[bin] != null) {
                    Iterator<ContiguousGenomicFeature> iteratorOverFeatures = setsOfFeaturesByGenomicBin[bin].iterator();
                    while (iteratorOverFeatures.hasNext()) {
                        ContiguousGenomicFeature feature = iteratorOverFeatures.next();
                        if (feature.getCoordinateOfEnd() < coordinate) {
                            if (setOfFeatures.size() == 0)
                                setOfFeatures.add(feature);
                            else {
                                int coordinateEndClosestFound = setOfFeatures.iterator().next().getCoordinateOfEnd();
                                if (feature.getCoordinateOfEnd() == coordinateEndClosestFound)
                                    setOfFeatures.add(feature);
                                else if (feature.getCoordinateOfEnd() > coordinateEndClosestFound) {
                                    setOfFeatures.clear();
                                    setOfFeatures.add(feature);
                                }
                            }
                        }
                    }
                }
                bin--;
            }
        }

        return setOfFeatures;
    }

    public HashSet<ContiguousGenomicFeature> getFirstFeatureWithStartCoordinateGreaterThan(String nameOfContig, int coordinate) {

        HashSet<ContiguousGenomicFeature> setOfFeatures = new HashSet<ContiguousGenomicFeature>();
        HashSet<ContiguousGenomicFeature> setsOfFeaturesByGenomicBin[] = this.mapContigNameToArrayOfBinnedFeatureSets.get(nameOfContig);
        if (setsOfFeaturesByGenomicBin != null) {
            int bin = Math.min(coordinate / this.sizeOfBins, setsOfFeaturesByGenomicBin.length -1);
            while (setOfFeatures.size() == 0 && bin < setsOfFeaturesByGenomicBin.length) {
                if (setsOfFeaturesByGenomicBin[bin] != null) {
                    Iterator<ContiguousGenomicFeature> iteratorOverFeatures = setsOfFeaturesByGenomicBin[bin].iterator();
                    while (iteratorOverFeatures.hasNext()) {
                        ContiguousGenomicFeature feature = iteratorOverFeatures.next();
                        if (feature.getCoordinateOfStart() > coordinate) {
                            if (setOfFeatures.size() == 0)
                                setOfFeatures.add(feature);
                            else {
                                int coordinateStartClosestFound = setOfFeatures.iterator().next().getCoordinateOfStart();
                                if (feature.getCoordinateOfStart() == coordinateStartClosestFound)
                                    setOfFeatures.add(feature);
                                else if (feature.getCoordinateOfStart() < coordinateStartClosestFound) {
                                    setOfFeatures.clear();
                                    setOfFeatures.add(feature);
                                }
                            }
                        }
                    }
                }
                bin++;
            }
        }

        return setOfFeatures;
    }

    public TreeSet<ContiguousGenomicFeature> getSetOfAllFeatures() {
        return this.setOfAllFeatures;
    }


    public static void main(String args[]) throws Exception  {

        ContiguousGenomicFeatureList listOfGenes = new ContiguousGenomicFeatureList(new File("/home/tuchbb/data/species/h_sapiens/ucsc_hg18_081115/gene_annotation.filtered.081115.gff"), 1000, 0, true, true, false, true);
        HashSet<ContiguousGenomicFeature> featuresOverlapping = listOfGenes.getSetOfFeaturesOverlapping("chr5", 132000000, 133000000).getFeaturesThatOverlapThisOneWithOverlapAnnotated();
        Iterator<ContiguousGenomicFeature> iterator = featuresOverlapping.iterator();
        while (iterator.hasNext()) {
            System.out.println(iterator.next());            
        }
    }

    /**
 * User: tuchbb
     * Date: Nov 25, 2008
     * Time: 9:28:22 AM
     * Revision: $Rev$
     * This code was originally developed as part of the SOLiD Whole Transcriptome package.
     */
    public static class OverlappingGenomicFeaturesSet {

        private int positionOfFeatureStart;
        private int positionOfFeatureEnd;

        private HashSet<ContiguousGenomicFeature> featuresThatOverlapThisOne;
        private HashSet<ContiguousGenomicFeature> featuresThatOverlapThisOneClonedWithOverlapFractions;

        private double maxOverlapFractionTheirPerspective;
        private double maxOverlapFractionMyPerspective;

        public OverlappingGenomicFeaturesSet(int positionOfFeatureStart, int positionOfFeatureEnd) {
            this.positionOfFeatureStart = positionOfFeatureStart;
            this.positionOfFeatureEnd = positionOfFeatureEnd;
            this.featuresThatOverlapThisOne = new HashSet<ContiguousGenomicFeature>();
            this.featuresThatOverlapThisOneClonedWithOverlapFractions = new HashSet<ContiguousGenomicFeature>();
            this.maxOverlapFractionMyPerspective = -1.0;
            this.maxOverlapFractionTheirPerspective = -1.0;
        }

        public void addOverlappingFeature(ContiguousGenomicFeature featureToAddWithOverlapAnnotated, ContiguousGenomicFeature featureToAddUniqueReference) {
            this.featuresThatOverlapThisOneClonedWithOverlapFractions.add(featureToAddWithOverlapAnnotated);
            this.featuresThatOverlapThisOne.add(featureToAddUniqueReference);
            this.maxOverlapFractionTheirPerspective = Math.max(this.maxOverlapFractionTheirPerspective, featureToAddWithOverlapAnnotated.getOverlapFractionFromItsPerspective());
            this.maxOverlapFractionMyPerspective = Math.max(this.maxOverlapFractionMyPerspective, featureToAddWithOverlapAnnotated.getOverlapFractionFromMyPerspective());
        }

        public int getPositionOfFeatureStart() {
            return positionOfFeatureStart;
        }

        public int getPositionOfFeatureEnd() {
            return positionOfFeatureEnd;
        }

        public HashSet<ContiguousGenomicFeature> getFeaturesThatOverlapThisOneWithOverlapAnnotated() {
            return featuresThatOverlapThisOneClonedWithOverlapFractions;
        }

        public HashSet<ContiguousGenomicFeature> getFeaturesThatOverlapThisOne() {
            return featuresThatOverlapThisOne;
        }

        public double getMaxOverlapFractionTheirPerspective() {
            return maxOverlapFractionTheirPerspective;
        }

        public double getMaxOverlapFractionMyPerspective() {
            return maxOverlapFractionMyPerspective;
        }

    }
}
