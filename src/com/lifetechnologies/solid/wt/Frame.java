package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Mar 3, 2009
 * Time: 12:10:08 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public enum Frame {

    NONE {
        public int toInt() { return -1; }
    },
    ZERO {
        public int toInt() { return 0; }
    },
    ONE {
        public int toInt() { return 1; }
    },
    TWO {
        public int toInt() { return 2; }
    };

    public abstract int toInt();

    public static Frame toFrame(int frame) throws Exception {
        if (frame == -1)   return NONE;
        else if (frame == 0)   return ZERO;
        else if (frame == 1)   return ONE;
        else if (frame == 2)   return TWO;
        else    throw new Exception("Frame int " + frame + " is unrecognized.");
    }


}
