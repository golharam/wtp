package com.lifetechnologies.solid.wt;

import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Dec 4, 2008
 * Time: 6:56:09 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 * 
 * Sort ContiguousGenomicFeatures by Reference ID, Start, End, Type, Strand
 */
public class ContiguousGenomicFeatureComparator implements Comparator<ContiguousGenomicFeature> {

	public static final ContiguousGenomicFeatureComparator INSTANCE = new ContiguousGenomicFeatureComparator();
	
    public int compare(ContiguousGenomicFeature contiguousGenomicFeatureA, ContiguousGenomicFeature contiguousGenomicFeatureB) {

        if (!contiguousGenomicFeatureA.equals(contiguousGenomicFeatureB)) {
            if (contiguousGenomicFeatureA.getIdOfReferenceSequence().equals(contiguousGenomicFeatureB.getIdOfReferenceSequence())) {

                if (contiguousGenomicFeatureA.getCoordinateOfStart() == contiguousGenomicFeatureB.getCoordinateOfStart()) {

                    if (contiguousGenomicFeatureA.getCoordinateOfEnd() == contiguousGenomicFeatureB.getCoordinateOfEnd()) {
                        if (contiguousGenomicFeatureA.getStrand() == contiguousGenomicFeatureA.getStrand()) {

                            if (contiguousGenomicFeatureA.getType() == null || contiguousGenomicFeatureB.getType() == null
                                    || contiguousGenomicFeatureA.getType().equalsIgnoreCase(contiguousGenomicFeatureB.getType())) {
                                return 0;
                            } else
                                return contiguousGenomicFeatureA.getType().compareTo(contiguousGenomicFeatureB.getType());
                            
                        } else if (contiguousGenomicFeatureA.getStrand() == Strand.POSITIVE )
                            return -1;
                        else
                            return 1;
                    } else if (contiguousGenomicFeatureA.getCoordinateOfEnd() < contiguousGenomicFeatureB.getCoordinateOfEnd())
                        return -1;
                    else
                        return 1;

                } else if (contiguousGenomicFeatureA.getCoordinateOfStart() < contiguousGenomicFeatureB.getCoordinateOfStart())
                    return -1;
                else
                    return 1;

            } else
                return contiguousGenomicFeatureA.getIdOfReferenceSequence().compareTo(contiguousGenomicFeatureB.getIdOfReferenceSequence());

        }

        return 0;
    }
}
