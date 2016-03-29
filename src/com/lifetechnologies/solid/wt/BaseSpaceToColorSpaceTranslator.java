package com.lifetechnologies.solid.wt;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * User: tuchbb
 * Date: Sep 24, 2008
 * Time: 12:47:40 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class BaseSpaceToColorSpaceTranslator {

    public static char[] NUCLEOTIDES_PLUS_IUPACS = {'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N'};

    private HashMap<String, String> mapTriNucleotideToDiColorsString;

    public BaseSpaceToColorSpaceTranslator() throws Exception {

        HashMap<String, HashSet<String>> mapTriNucleotideToDiColorSet = new HashMap<String, HashSet<String>>();

        for (int i = 0; i < NUCLEOTIDES_PLUS_IUPACS.length; i++)
            for (int j = 0; j < NUCLEOTIDES_PLUS_IUPACS.length; j++)
                for (int k = 0; k < NUCLEOTIDES_PLUS_IUPACS.length; k++) {
                    String triNucleotide = NUCLEOTIDES_PLUS_IUPACS[i] + "" + NUCLEOTIDES_PLUS_IUPACS[j] + "" + NUCLEOTIDES_PLUS_IUPACS[k];
                    mapTriNucleotideToDiColorSet.put(triNucleotide, new HashSet<String>());
                    char basesAtFirstPosition[] = expandIUPAC(NUCLEOTIDES_PLUS_IUPACS[i]);
                    char basesAtSecondPosition[] = expandIUPAC(NUCLEOTIDES_PLUS_IUPACS[j]);
                    char basesAtThirdPosition[] = expandIUPAC(NUCLEOTIDES_PLUS_IUPACS[k]);
                    for (int l = 0; l < basesAtFirstPosition.length; l++)
                        for (int m = 0; m < basesAtSecondPosition.length; m++)
                            for (int n = 0; n < basesAtThirdPosition.length; n++)
                                mapTriNucleotideToDiColorSet.get(triNucleotide).add(translatePairOfBases(basesAtFirstPosition[l], basesAtSecondPosition[m])
                                                                             + "" + translatePairOfBases(basesAtSecondPosition[m], basesAtThirdPosition[n]));

                }

        this.mapTriNucleotideToDiColorsString = new HashMap<String, String>();
        Iterator<String> iteratorTriNucleotides = mapTriNucleotideToDiColorSet.keySet().iterator();
        while (iteratorTriNucleotides.hasNext()) {
            String triNucleotide = iteratorTriNucleotides.next();
            String diColorString = "";
            Iterator<String> iteratorDiColors = mapTriNucleotideToDiColorSet.get(triNucleotide).iterator();
            while (iteratorDiColors.hasNext()) {
                diColorString += iteratorDiColors.next() + "|";
            }
            this.mapTriNucleotideToDiColorsString.put(triNucleotide, diColorString.substring(0, diColorString.length() -1));            
        }

    }

    public String[] translateSequenceWithIUPACs(String sequenceInBaseSpace) throws Exception {

        String sequenceInColorSpace[] = new String[sequenceInBaseSpace.length() -2];

        char sequenceInBaseSpaceAsCharArray[] = sequenceInBaseSpace.toCharArray();
        for (int i = 0; i < sequenceInColorSpace.length; i++)
            sequenceInColorSpace[i] = this.mapTriNucleotideToDiColorsString.get(sequenceInBaseSpaceAsCharArray[i] + "" + sequenceInBaseSpaceAsCharArray[i +1] + "" + sequenceInBaseSpaceAsCharArray[i +2]);

        return sequenceInColorSpace;
    }

    private static char[] expandIUPAC(char iupac) {
        if (iupac == 'A')
            return new char[]{'A'};
        else if (iupac == 'C')
            return new char[]{'C'};
        else if (iupac == 'G')
            return new char[]{'G'};
        else if (iupac == 'T')
            return new char[]{'T'};
        else if (iupac == 'N')
            return new char[]{'A', 'C', 'G', 'T'};
        else if (iupac == 'R')
            return new char[]{'A', 'G'};
        else if (iupac == 'Y')
            return new char[]{'C', 'T'};
        else if (iupac == 'M')
            return new char[]{'A', 'C'};
        else if (iupac == 'K')
            return new char[]{'G', 'T'};
        else if (iupac == 'W')
            return new char[]{'A', 'T'};
        else if (iupac == 'S')
            return new char[]{'C', 'G'};
        else if (iupac == 'B')
            return new char[]{'C', 'G', 'T'};
        else if (iupac == 'D')
            return new char[]{'A', 'G', 'T'};
        else if (iupac == 'H')
            return new char[]{'A', 'C', 'T'};
        else if (iupac == 'V')
            return new char[]{'A', 'C', 'G'};
        else
            return null;

    }

    public static String translateSequence(String sequenceInBaseSpace) throws Exception {

        String sequenceInColorSpace = "";

        char sequenceInBaseSpaceAsCharArray[] = sequenceInBaseSpace.toCharArray();
        for (int i = 0; i < sequenceInBaseSpaceAsCharArray.length -1; i++)
            sequenceInColorSpace += translatePairOfBases(sequenceInBaseSpaceAsCharArray[i], sequenceInBaseSpaceAsCharArray[i +1]);

        return sequenceInColorSpace;
    }

    private static int translatePairOfBases(char baseFivePrime, char baseThreePrime) throws Exception {

        String pairOfBases = "" + baseFivePrime + baseThreePrime;   //).toUpperCase();
        if (pairOfBases.matches("AA|CC|GG|TT"))
            return 0;
        else if (pairOfBases.matches("AC|CA|GT|TG"))
            return 1;
        else if (pairOfBases.matches("AG|GA|CT|TC"))
            return 2;
        else if (pairOfBases.matches("AT|TA|CG|GC"))
            return 3;
        else
            throw new Exception("Pair of bases is not valid: " + pairOfBases);

    }

    public static void main(String args[]) throws Exception {
        new BaseSpaceToColorSpaceTranslator();
        System.out.println(translatePairOfBases('T', 'A'));
        double temp = 1e-100;
        System.out.println(temp);
        
    }
}
