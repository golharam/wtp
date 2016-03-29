package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Dec 15, 2008
 * Time: 8:50:44 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class SequenceHeaderPair {
    private String header;
    private String sequence;

    public SequenceHeaderPair(String header, String sequence) {
        this.header = header;
        this.sequence = sequence;
    }

    public String getHeader() {
        return header;
    }

    public void setHeader(String header) {
        this.header = header;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
}
