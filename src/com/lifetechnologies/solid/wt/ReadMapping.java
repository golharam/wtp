package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Aug 20, 2008
 * Time: 3:32:59 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ReadMapping {

    private String idOfRead;
    private int indexOfMatchingReferenceSequence;
    private String idOfMatchingReferenceSequence;
    private int positionOfAlignmentStartInReferenceSequence;
    private int numberOfMismatches;

    public ReadMapping(String idOfRead, String idOfMatchingReferenceSequence, int indexOfMatchingReferenceSequence, int positionOfAlignmentStartInReferenceSequence, int numberOfMismatches) {
        this.idOfRead = idOfRead;
        this.idOfMatchingReferenceSequence = idOfMatchingReferenceSequence;
        this.indexOfMatchingReferenceSequence = indexOfMatchingReferenceSequence;
        this.positionOfAlignmentStartInReferenceSequence = positionOfAlignmentStartInReferenceSequence;
        this.numberOfMismatches = numberOfMismatches;
    }

    public ReadMapping(String idOfMatchingReferenceSequence, int positionOfAlignmentStartInReferenceSequence, int numberOfMismatches) {
        this.idOfMatchingReferenceSequence = idOfMatchingReferenceSequence;
        this.positionOfAlignmentStartInReferenceSequence = positionOfAlignmentStartInReferenceSequence;
        this.numberOfMismatches = numberOfMismatches;
    }

    public ReadMapping(int indexOfMatchingReferenceSequence, int positionOfAlignmentStartInReferenceSequence) {
        this.indexOfMatchingReferenceSequence = indexOfMatchingReferenceSequence;
        this.positionOfAlignmentStartInReferenceSequence = positionOfAlignmentStartInReferenceSequence;
    }

    public int getIndexOfMatchingReferenceSequence() {
        return indexOfMatchingReferenceSequence;
    }

    public void setIndexOfMatchingReferenceSequence(int indexOfMatchingReferenceSequence) {
        this.indexOfMatchingReferenceSequence = indexOfMatchingReferenceSequence;
    }

    public String getIdOfMatchingReferenceSequence() {
        return idOfMatchingReferenceSequence;
    }

    public void setIdOfMatchingReferenceSequence(String idOfMatchingReferenceSequence) {
        this.idOfMatchingReferenceSequence = idOfMatchingReferenceSequence;
    }

    public int getPositionOfAlignmentStartInReferenceSequence() {
        return positionOfAlignmentStartInReferenceSequence;
    }

    public void setPositionOfAlignmentStartInReferenceSequence(int positionOfAlignmentStartInReferenceSequence) {
        this.positionOfAlignmentStartInReferenceSequence = positionOfAlignmentStartInReferenceSequence;
    }

    public int getNumberOfMismatches() {
        return numberOfMismatches;
    }

    public void setNumberOfMismatches(int numberOfMismatches) {
        this.numberOfMismatches = numberOfMismatches;
    }

    public String getIdOfRead() {
        return idOfRead;
    }

    public void setIdOfRead(String idOfRead) {
        this.idOfRead = idOfRead;
    }


}
