package com.lifetechnologies.solid.wt;

import com.lifetechnologies.solid.wt.Utilities;

import java.io.File;
import java.io.IOException;
import java.util.TreeMap;

/**
 * User: btuch
 * Date: Feb 16, 2005
 */
public class CodonTable {

    private TreeMap<String, String> mapCodonToAminoAcid = new TreeMap<String, String>();

    public CodonTable(File codonTableFile) throws IOException {
        this.mapCodonToAminoAcid = Utilities.loadStringToStringMapFromFile(codonTableFile, 0, 2, "\t", 0);
    }

    public String translateCodon(String codon) {

        String codonUpperCase = codon.toUpperCase();
        if (codonUpperCase.indexOf('N') > -1)   return "X";

        return this.mapCodonToAminoAcid.get(codonUpperCase);
    }

    public String translateSequence(String sequenceORF) {
        String seqeunceTranslatedORF = "";
        for (int i = 0; i < sequenceORF.length(); i=i+3)
            seqeunceTranslatedORF += translateCodon(sequenceORF.substring(i, i+3));

        return seqeunceTranslatedORF.toUpperCase();
    }
}
