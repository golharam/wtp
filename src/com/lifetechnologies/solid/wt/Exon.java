package com.lifetechnologies.solid.wt;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.Strand;

/**
 * User: tuchbb
 * Date: Mar 3, 2009
 * Time: 9:03:20 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class Exon extends ContiguousGenomicFeature {

    private int positionOfCodingSequenceStart;
    private int positionOfCodingSequenceEnd;
    private Frame frameOfCodingSequence;

    private Exon previousExon;
    private Exon nextExon;

    public Exon(String idOfReferenceSequence, Strand strand, int coordinateOfStart, int coordinateOfEnd) {
        super(idOfReferenceSequence, strand, coordinateOfStart, coordinateOfEnd);
    }

    public Exon getPreviousExon() {
        return previousExon;
    }

    public void setPreviousExon(Exon previousExon) {
        this.previousExon = previousExon;
    }

    public Exon getNextExon() {
        return nextExon;
    }

    public void setNextExon(Exon nextExon) {
        this.nextExon = nextExon;
    }

    public int getPositionOfCodingSequenceStart() {
        return positionOfCodingSequenceStart;
    }

    public void setPositionOfCodingSequenceStart(int positionOfCodingSequenceStart) {
        this.positionOfCodingSequenceStart = positionOfCodingSequenceStart;
    }

    public int getPositionOfCodingSequenceEnd() {
        return positionOfCodingSequenceEnd;
    }

    public void setPositionOfCodingSequenceEnd(int positionOfCodingSequenceEnd) {
        this.positionOfCodingSequenceEnd = positionOfCodingSequenceEnd;
    }

    public Frame getFrameOfCodingSequence() {
        return frameOfCodingSequence;
    }

    public void setFrameOfCodingSequence(Frame frameOfCodingSequence) {
        this.frameOfCodingSequence = frameOfCodingSequence;
    }

}
