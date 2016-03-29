package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Oct 17, 2008
 * Time: 9:56:28 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class Read {

    private String nameOfRead;
    private String sequenceOfRead;

    public Read(String nameOfRead, String sequenceOfRead) {
        this.nameOfRead = nameOfRead;
        this.sequenceOfRead = sequenceOfRead;
    }

    public String getNameOfRead() {
        return nameOfRead;
    }

    public void setNameOfRead(String nameOfRead) {
        this.nameOfRead = nameOfRead;
    }

    public String getSequenceOfRead() {
        return sequenceOfRead;
    }

    public void setSequenceOfRead(String sequenceOfRead) {
        this.sequenceOfRead = sequenceOfRead;
    }
}
