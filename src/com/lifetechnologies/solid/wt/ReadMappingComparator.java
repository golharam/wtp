package com.lifetechnologies.solid.wt;

import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Oct 17, 2008
 * Time: 9:38:44 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */

public class ReadMappingComparator implements Comparator<ReadMapping> {


    public int compare(ReadMapping readMappingA, ReadMapping readMappingB) {

        if (!readMappingA.equals(readMappingB)) {
            if (readMappingA.getIndexOfMatchingReferenceSequence() == readMappingB.getIndexOfMatchingReferenceSequence()) {

                int positionAbsoluteA = Math.abs(readMappingA.getPositionOfAlignmentStartInReferenceSequence());
                int positionAbsoluteB = Math.abs(readMappingB.getPositionOfAlignmentStartInReferenceSequence());
                if (positionAbsoluteA == positionAbsoluteB)
                    return readMappingA.getIdOfRead().compareTo(readMappingB.getIdOfRead());
                else if (positionAbsoluteA > positionAbsoluteB)
                    return 1;
                else
                    return -1;

            } else if (readMappingA.getIndexOfMatchingReferenceSequence() > readMappingB.getIndexOfMatchingReferenceSequence())
                return 1;
            else
                return -1;
        }

        return 0;
    }
}