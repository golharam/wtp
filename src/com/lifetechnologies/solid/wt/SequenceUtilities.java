package com.lifetechnologies.solid.wt;

/**
 * User: tuchbb
 * Date: Mar 5, 2009
 * Time: 9:25:21 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class SequenceUtilities {

    public static String reverseComplement(String sequence) {

        String reverseComplement = "";
        for (int i = sequence.length() - 1; i >= 0; i--)
            reverseComplement += complement(sequence.charAt(i));

        return reverseComplement;
    }

    public static StringBuffer reverseComplement(StringBuffer sequence) {

        StringBuffer reverseComplement = new StringBuffer();
        for (int i = sequence.length() - 1; i >= 0; i--)
            reverseComplement.append(complement(sequence.charAt(i)));

        return reverseComplement;
    }

    public static StringBuffer complement(StringBuffer sequence) {

        StringBuffer complement = new StringBuffer();
        for (int i = 0; i < sequence.length(); i++)
            complement.append(complement(sequence.charAt(i)));

        return complement;
    }

    private static char complement(char c) {
        char base = Character.toLowerCase(c);
        if (base == 'a')
            return 't';
        else if (base == 't')
            return 'a';
        else if (base == 'g')
            return 'c';
        else if (base == 'c')
            return 'g';
        else if (base == 'm')
            return 'k';
        else if (base == 'r')
            return 'y';
        else if (base == 'w')
            return 'w';
        else if (base == 's')
            return 's';
        else if (base == 'y')
            return 'r';
        else if (base == 'k')
            return 'm';
        else if (base == 'v')
            return 'b';
        else if (base == 'h')
            return 'd';
        else if (base == 'd')
            return 'h';
        else if (base == 'b')
            return 'v';
        else if (base == '-')   // added 072406
            return '-';
        else
            return 'n';
    }

}
