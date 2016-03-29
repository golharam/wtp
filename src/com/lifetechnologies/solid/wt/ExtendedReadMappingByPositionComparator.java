package com.lifetechnologies.solid.wt;

import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Dec 15, 2008
 * Time: 2:10:08 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExtendedReadMappingByPositionComparator implements Comparator<ExtendedReadMapping> {

	public static final Comparator<ExtendedReadMapping> INSTANCE = new ExtendedReadMappingByPositionComparator();

    public int compare(ExtendedReadMapping extendedReadMappingA, ExtendedReadMapping extendedReadMappingB) {

    	if (extendedReadMappingA == null) {
    		if (extendedReadMappingB == null) return 0;
    		return -1;
    	}
    	if (extendedReadMappingB == null) return 1;
    	
        if (!extendedReadMappingA.equals(extendedReadMappingB)) {
            if (extendedReadMappingA.getIndexOfMatchingReferenceSequence() == extendedReadMappingB.getIndexOfMatchingReferenceSequence()) {

                int positionAbsoluteA = Math.abs(extendedReadMappingA.getPositionOfAlignmentStartInReferenceSequence());
                int positionAbsoluteB = Math.abs(extendedReadMappingB.getPositionOfAlignmentStartInReferenceSequence());

                if (positionAbsoluteA == positionAbsoluteB)
                    return extendedReadMappingA.getIdOfRead().compareTo(extendedReadMappingB.getIdOfRead());
                else if (positionAbsoluteA > positionAbsoluteB)
                    return 1;
                else
                    return -1;

            } else if (extendedReadMappingA.getIndexOfMatchingReferenceSequence() > extendedReadMappingB.getIndexOfMatchingReferenceSequence())
                return 1;
            else
                return -1;
        }

        return 0;
    }

}
