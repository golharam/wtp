package com.lifetechnologies.solid.wt;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.Strand;

import java.io.*;

/**
 * User: tuchbb
 * Date: Dec 22, 2008
 * Time: 2:21:27 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class SNPTableReader {

    File fileSortedSNPTable;

    private BufferedReader readerSortedSNPTable; // = new BufferedReader(new FileReader(args[2]));
    private ContiguousGenomicFeature featureCurrent;
    private Character nucleotideInNCBIReference;
    private Character nucleotideInUCSCReference;
    private String nucleotidesObserved[];
    private boolean peeked = false;

    public SNPTableReader(File fileSortedSNPTable) throws FileNotFoundException {
        this.fileSortedSNPTable = fileSortedSNPTable;
        this.readerSortedSNPTable = new BufferedReader(new FileReader(fileSortedSNPTable));
    }

    protected void finalize() throws Throwable {
        if (this.readerSortedSNPTable != null)
            this.readerSortedSNPTable.close();
        super.finalize();
    }

    public void reset() throws IOException {
        this.readerSortedSNPTable = new BufferedReader(new FileReader(this.fileSortedSNPTable));
    }

    public ContiguousGenomicFeature peekAtNextEntry() throws Exception {
        if (!this.peeked) {
            String line = this.readerSortedSNPTable.readLine();
            while (line != null && line.startsWith("#"))
                line = this.readerSortedSNPTable.readLine();

            if (line != null) {

                String tokens[] = line.split("\t");
                this.featureCurrent = new ContiguousGenomicFeature(tokens[1],  Strand.toStrand(tokens[6].toCharArray()[0]), Integer.parseInt(tokens[2]), Integer.parseInt(tokens[3]));

                this.featureCurrent.setLabel(tokens[4]);    // SNP name
                this.featureCurrent.setType(tokens[15]);    // functional category
                this.nucleotideInNCBIReference = tokens[7].toCharArray()[0];   // NCBI reference nucleotide
                this.nucleotideInUCSCReference = tokens[8].toCharArray()[0];   // UCSC reference nucleotide
                this.nucleotidesObserved = tokens[9].split("/");   // observed nucleotides
                if (this.featureCurrent.getStrand() == Strand.NEGATIVE)
                    for (int i = 0; i < nucleotidesObserved.length; i++)
                        this.nucleotidesObserved[i] = complement(this.nucleotidesObserved[i]);


            } else {
                this.featureCurrent = null;
                this.nucleotideInNCBIReference = null;
                this.nucleotideInUCSCReference = null;
                this.nucleotidesObserved = null;
            }
            this.peeked = true;
        }
        return this.featureCurrent;
    }

    private String complement(String nucleotide) {
        if (nucleotide.equalsIgnoreCase("A"))
            return "T";
        if (nucleotide.equalsIgnoreCase("C"))
            return "G";
        if (nucleotide.equalsIgnoreCase("G"))
            return "C";
        if (nucleotide.equalsIgnoreCase("T"))
            return "A";
        return null;
    }

    public ContiguousGenomicFeature nextEntry() throws Exception {
        if (!this.peeked)
            peekAtNextEntry();
        this.peeked = false;
        return this.featureCurrent;
    }

    public Character getNucleotideInNCBIReference() {
        return nucleotideInNCBIReference;
    }

    public Character getNucleotideInUCSCReference() {
        return nucleotideInUCSCReference;
    }

    public String[] getNucleotidesObserved() {
        return nucleotidesObserved;
    }

}
