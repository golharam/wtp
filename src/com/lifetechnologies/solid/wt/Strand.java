package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
* Date: Feb 13, 2009
* Time: 8:21:53 AM
* Revision: $Rev$
* This code was originally developed as part of the SOLiD Whole Transcriptome package.
*/
public enum Strand {

    POSITIVE {
        public char toChar() { return '+'; }
    } ,
    NEGATIVE {
        public char toChar() { return '-'; }
    },
    EITHER {
        public char toChar() { return '.'; }
    };

    public abstract char toChar();
    //abstract boolean isMatch(Strand strand);

    public static Strand toStrand(char character) throws IllegalArgumentException {
        if (character == '+')   return POSITIVE;
        else if (character == '-')   return NEGATIVE;
        else if (character == '.')   return EITHER;
        else    throw new IllegalArgumentException("Strand character " + character + " is unrecognized.");
    }
}
