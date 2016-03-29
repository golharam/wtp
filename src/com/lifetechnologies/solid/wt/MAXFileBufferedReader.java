package com.lifetechnologies.solid.wt;

import java.io.*;

/**
 * User: tuchbb
 * Date: Nov 19, 2008
 * Time: 3:39:21 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class MAXFileBufferedReader {

    private BufferedReader readerMAXFile; // = new BufferedReader(new FileReader(args[2]));
    private ExtendedReadMapping[] extendedReadMappings;
    private String idOfCurrentRead;
    private String sequenceOfCurrentRead;
    private boolean peeked = false;

    public MAXFileBufferedReader(File fileMAX) throws FileNotFoundException {
        this.readerMAXFile = new BufferedReader(new FileReader(fileMAX));
    }

    protected void finalize() throws Throwable {
        if (this.readerMAXFile != null)
            this.readerMAXFile.close();
        super.finalize();
    }

    public ExtendedReadMapping[] peekAtNextEntry() throws IOException {
        if (!this.peeked) {
            String header = this.readerMAXFile.readLine();
            if (header != null) {
                this.extendedReadMappings = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(header.substring(1));
                int indexOfCommaInHeader = header.indexOf(',');
                if (indexOfCommaInHeader > -1)
                    this.idOfCurrentRead = header.substring(1, indexOfCommaInHeader);
                else
                    this.idOfCurrentRead = header.substring(1);
                this.sequenceOfCurrentRead = this.readerMAXFile.readLine();
            } else {
                this.extendedReadMappings = null;
                this.idOfCurrentRead = null;
                this.sequenceOfCurrentRead = null;
            }
            this.peeked = true;
        }
        return this.extendedReadMappings;
    }

    public ExtendedReadMapping[] nextEntry() throws IOException {
        if (!this.peeked)
            peekAtNextEntry();
        this.peeked = false;
        return this.extendedReadMappings;
    }

    public String getIdOfCurrentRead() {
        return idOfCurrentRead;
    }

    public String getSequenceOfCurrentRead() {
        return sequenceOfCurrentRead;
    }
}
