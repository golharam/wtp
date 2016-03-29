package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Oct 1, 2008
 * Time: 11:39:44 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class Insertion {

    private int position;
    private int length;

    public Insertion(int position, int length) {
        this.position = position;
        this.length = length;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public String toString() {
        return position + "." + length;
    }
}
