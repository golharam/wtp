package com.lifetechnologies.solid.wt;

import java.util.HashMap;
import java.util.HashSet;

/**
 * User: tuchbb
 * Date: Nov 25, 2008
 * Time: 9:28:22 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class OverlappingGenomicFeaturesSet {

    private int positionOfFeatureStart;
    private int positionOfFeatureEnd;

    private HashSet<ContiguousGenomicFeature> featuresThatOverlapThisOne;

    private double maxOverlapFractionTheirPerspective;
    private double maxOverlapFractionMyPerspective;

    public OverlappingGenomicFeaturesSet(int positionOfFeatureStart, int positionOfFeatureEnd) {
        this.positionOfFeatureStart = positionOfFeatureStart;
        this.positionOfFeatureEnd = positionOfFeatureEnd;
        this.featuresThatOverlapThisOne = new HashSet<ContiguousGenomicFeature>();
        this.maxOverlapFractionMyPerspective = -1.0;
        this.maxOverlapFractionTheirPerspective = -1.0;
    }

    public void addOverlappingFeature(ContiguousGenomicFeature featureToAdd) {
        this.featuresThatOverlapThisOne.add(featureToAdd);
        this.maxOverlapFractionTheirPerspective = Math.max(this.maxOverlapFractionTheirPerspective, featureToAdd.getOverlapFractionFromItsPerspective());
        this.maxOverlapFractionMyPerspective = Math.max(this.maxOverlapFractionMyPerspective, featureToAdd.getOverlapFractionFromMyPerspective());        
    }

    public int getPositionOfFeatureStart() {
        return positionOfFeatureStart;
    }

    public int getPositionOfFeatureEnd() {
        return positionOfFeatureEnd;
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

    public static OverlappingGenomicFeaturesSet findOverlappingGenomicFeatures(HashMap<Integer, ContiguousGenomicFeature[]> mapPositionToFeatures,
                                                                               ContiguousGenomicFeature featureToExamine, double minOverlapFraction) {

        OverlappingGenomicFeaturesSet featuresThatOverlap = new OverlappingGenomicFeaturesSet(featureToExamine.getCoordinateOfStart(), featureToExamine.getCoordinateOfEnd());

        HashSet<ContiguousGenomicFeature> setOfFeatureAlreadyExamined = new HashSet<ContiguousGenomicFeature>();
        for (int positionCurrent = featureToExamine.getCoordinateOfStart(); positionCurrent <= featureToExamine.getCoordinateOfEnd(); positionCurrent++) {

            if (mapPositionToFeatures.containsKey(positionCurrent)) {
                ContiguousGenomicFeature featuresOverlapingCurrentPosition[] = mapPositionToFeatures.get(positionCurrent);
                for (int i = 0; i < featuresOverlapingCurrentPosition.length; i++)
                    if (!setOfFeatureAlreadyExamined.contains(featuresOverlapingCurrentPosition[i])) {
                        setOfFeatureAlreadyExamined.add(featuresOverlapingCurrentPosition[i]);
                        int lengthOfOverlapingSequence = 1 + Math.min(featureToExamine.getCoordinateOfEnd(), featuresOverlapingCurrentPosition[i].getCoordinateOfEnd())
                                                           - Math.max(featureToExamine.getCoordinateOfStart(), featuresOverlapingCurrentPosition[i].getCoordinateOfStart());
                        double overlapFromPredictedFeaturePerspective = lengthOfOverlapingSequence / (double)featureToExamine.getLength();
                        double overlapFromAnnotatedFeaturePerspective = lengthOfOverlapingSequence / (double)featuresOverlapingCurrentPosition[i].getLength();

                        ContiguousGenomicFeature featureAnnotatedWithOverlap = new ContiguousGenomicFeature(featuresOverlapingCurrentPosition[i].getCoordinateOfStart(), featuresOverlapingCurrentPosition[i].getCoordinateOfEnd());
                        featureAnnotatedWithOverlap.setOverlapFractionFromMyPerspective(overlapFromAnnotatedFeaturePerspective);
                        featureAnnotatedWithOverlap.setOverlapFractionFromItsPerspective(overlapFromPredictedFeaturePerspective);
                        featureAnnotatedWithOverlap.setParent(featuresOverlapingCurrentPosition[i]);
                        if (featuresOverlapingCurrentPosition[i].getLabel() != null)
                            featureAnnotatedWithOverlap.setLabel(featuresOverlapingCurrentPosition[i].getLabel());
                        if (overlapFromPredictedFeaturePerspective >= minOverlapFraction || overlapFromAnnotatedFeaturePerspective >= minOverlapFraction)
                            featuresThatOverlap.addOverlappingFeature(featureAnnotatedWithOverlap);

                    }
            }
        }

        /*System.out.println("mapPositionToFeatures.size():\t" + mapPositionToFeatures.size());
        System.out.println("featureToExamine:\t" + featureToExamine);
        System.out.println("minOverlapFraction:\t" + minOverlapFraction);
        System.out.println("featuresThatOverlap.getFeaturesThatOverlapThisOne().size():\t" + featuresThatOverlap.getFeaturesThatOverlapThisOne().size());
          */
        return featuresThatOverlap;
    }
}
