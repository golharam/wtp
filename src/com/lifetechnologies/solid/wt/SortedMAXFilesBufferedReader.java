package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.TreeMap;
import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Nov 28, 2008
 * Time: 11:55:31 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class SortedMAXFilesBufferedReader {

    private static String PREFIX_MAX_FILE_INDEX = "file_";

    private File filesMAXSorted[];
    private Comparator<ExtendedReadMapping> comparatorOfMAXEntries;
    private BufferedReader readersMAXSorted[];
    private TreeMap<ExtendedReadMapping, Integer> mapExtendedReadMappingToMAXFileIndex;

    private boolean storeReadSequence;

    private long numberOfReadsProcessed;

    public SortedMAXFilesBufferedReader(File[] filesMAXSorted, Comparator<ExtendedReadMapping> comparatorOfMAXEntries, boolean storeReadSequence) throws IOException {

        this.storeReadSequence = storeReadSequence;
        this.filesMAXSorted = filesMAXSorted;
        this.comparatorOfMAXEntries = comparatorOfMAXEntries;

        initialize();
    }

    private void initialize() throws IOException {

        this.numberOfReadsProcessed = 0;

        this.readersMAXSorted = new BufferedReader[this.filesMAXSorted.length];
        for (int i = 0; i < this.filesMAXSorted.length; i++)
            this.readersMAXSorted[i] = new BufferedReader(new FileReader(this.filesMAXSorted[i]));

        this.mapExtendedReadMappingToMAXFileIndex = new TreeMap<ExtendedReadMapping, Integer>(this.comparatorOfMAXEntries);
        for (int indexMAXFile = 0; indexMAXFile < this.readersMAXSorted.length; indexMAXFile++) {
            String header = this.readersMAXSorted[indexMAXFile].readLine();
            if (header != null) {
                ExtendedReadMapping extendedReadMappings[] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(header.substring(1));
                extendedReadMappings[0].setIdOfRead(PREFIX_MAX_FILE_INDEX + indexMAXFile + "_" + extendedReadMappings[0].getIdOfRead());
                String sequenceOfRead = this.readersMAXSorted[indexMAXFile].readLine();
                if (this.storeReadSequence)
                    extendedReadMappings[0].setSequenceOfRead(sequenceOfRead);
                this.mapExtendedReadMappingToMAXFileIndex.put(extendedReadMappings[0], indexMAXFile);
            }
        }
    }

    protected void finalize() throws Throwable {

        for (int i = 0; i < this.readersMAXSorted.length; i++) {
            if (this.readersMAXSorted[i] != null)
                this.readersMAXSorted[i].close();
        }
        super.finalize();
    }



    public ExtendedReadMapping peekAtNextEntry() throws IOException {
        if (this.mapExtendedReadMappingToMAXFileIndex.size() > 0 )
            return this.mapExtendedReadMappingToMAXFileIndex.firstKey();
        return null;
    }

    public ExtendedReadMapping nextEntry() throws IOException {

        if (this.mapExtendedReadMappingToMAXFileIndex.size() > 0) {

            ExtendedReadMapping extendedReadMappingFirst = this.mapExtendedReadMappingToMAXFileIndex.firstKey();
            Integer indexOfCorrespondingPartition = this.mapExtendedReadMappingToMAXFileIndex.get(extendedReadMappingFirst);
            this.mapExtendedReadMappingToMAXFileIndex.remove(extendedReadMappingFirst);

            String header = this.readersMAXSorted[indexOfCorrespondingPartition].readLine();
            if (header != null) {
                ExtendedReadMapping extendedReadMappings[] = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(header.substring(1));
                extendedReadMappings[0].setIdOfRead(PREFIX_MAX_FILE_INDEX + indexOfCorrespondingPartition + "_" + extendedReadMappings[0].getIdOfRead());
                String sequenceOfRead = this.readersMAXSorted[indexOfCorrespondingPartition].readLine();
                if (this.storeReadSequence)
                    extendedReadMappings[0].setSequenceOfRead(sequenceOfRead);
                this.mapExtendedReadMappingToMAXFileIndex.put(extendedReadMappings[0], indexOfCorrespondingPartition);
            }

            this.numberOfReadsProcessed++;

            extendedReadMappingFirst.setIdOfRead(extendedReadMappingFirst.getIdOfRead().substring(extendedReadMappingFirst.getIdOfRead().indexOf('_', PREFIX_MAX_FILE_INDEX.length())) +1);

            return extendedReadMappingFirst;

        }

        return null;
    }

    public long getNumberOfReadsProcessed() {
        return this.numberOfReadsProcessed;
    }

    public File[] getFilesMAXSorted() {
        return filesMAXSorted;
    }

    public void reset() throws IOException {
        initialize();        
    }
}
