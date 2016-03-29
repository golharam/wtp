package com.lifetechnologies.solid.wt;

import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Oct 2, 2008
 * Time: 10:33:37 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExtendedReadMappingByScoreComparator implements Comparator<ExtendedReadMapping> {

	public static final ExtendedReadMappingByScoreComparator INSTANCE = new ExtendedReadMappingByScoreComparator();
	
    public int compare(ExtendedReadMapping extendedReadMappingA, ExtendedReadMapping extendedReadMappingB) {

        if (!extendedReadMappingA.equals(extendedReadMappingB)) {
            if (extendedReadMappingA.getScore() == extendedReadMappingB.getScore()) {

                if (extendedReadMappingA.getIndexOfMatchingReferenceSequence() == extendedReadMappingB.getIndexOfMatchingReferenceSequence()) {
                    int positionAbsoluteA = Math.abs(extendedReadMappingA.getPositionOfAlignmentStartInReferenceSequence());
                    int positionAbsoluteB = Math.abs(extendedReadMappingB.getPositionOfAlignmentStartInReferenceSequence());
                    return positionAbsoluteA - positionAbsoluteB;
                } else if (extendedReadMappingA.getIndexOfMatchingReferenceSequence() > extendedReadMappingB.getIndexOfMatchingReferenceSequence())
                    return 1;
                else
                    return -1;

            } else if (extendedReadMappingA.getScore() > extendedReadMappingB.getScore())
                return -1;
            else
                return 1;
        }

        return 0;
    }
}
