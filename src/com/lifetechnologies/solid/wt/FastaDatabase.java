package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

import com.lifetechnologies.util.FileHeader;

/**
 * User: tuchbb
 * Date: Aug 13, 2008
 * Time: 10:53:19 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class FastaDatabase {

    protected File fileFastaDatabase = null;
    private int numberOfSequencesInDB;
    private ArrayList<String> listOfHeaders = null;
    protected HashMap<String, String> mapOfHeadersToSequences = null;
    private BufferedReader readerFastaFile = null;
    private String headerCurrentSequence;

    private static Logger logger = Logger.getLogger(FastaDatabase.class.getSimpleName());
    
    protected FastaDatabase() {}
    
    public FastaDatabase(File fastaDatabaseFile) throws IOException {
        this.fileFastaDatabase = fastaDatabaseFile;
        initialize();
    }

    private void initialize() throws IOException {
        this.numberOfSequencesInDB = -1;
        this.readerFastaFile = new BufferedReader(new FileReader(this.fileFastaDatabase));
        this.headerCurrentSequence = "";
        while (this.headerCurrentSequence != null && !this.headerCurrentSequence.startsWith(">"))
            this.headerCurrentSequence = this.readerFastaFile.readLine();
    }

    protected void finalize() throws Throwable {
        if (this.readerFastaFile != null)
            this.readerFastaFile.close();
        super.finalize();
    }

    public void reset() throws IOException {
        initialize();
    }
    
    public SequenceHeaderPair nextSequence(boolean toUpperCase) throws Exception {

        SequenceHeaderPair sequenceHeaderPair = null;
        if (this.headerCurrentSequence != null) {

            StringBuffer stringBufferSequence = new StringBuffer();
            String line;
            while ((line = this.readerFastaFile.readLine()) != null && !line.startsWith(">"))
                stringBufferSequence.append(line.trim());
            sequenceHeaderPair = new SequenceHeaderPair(this.headerCurrentSequence, stringBufferSequence.toString());
            if (toUpperCase)    sequenceHeaderPair.setSequence(sequenceHeaderPair.getSequence().toUpperCase());
            this.headerCurrentSequence = line;
        }

        return sequenceHeaderPair;
    }

    public int getNumberOfSequencesInDatabase() throws Exception {
        if (this.numberOfSequencesInDB == -1) {
            this.numberOfSequencesInDB = 0;
            BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
            String line = null;
            while ((line = reader.readLine()) != null)
                if (line.startsWith(">"))  this.numberOfSequencesInDB++;
            reader.close();
        }
        return numberOfSequencesInDB;
    }

    @SuppressWarnings("unchecked")
    public List<String> getListOfHeaders(boolean storeCopy, boolean truncateAfterFirstSpace) throws Exception {
        if (this.listOfHeaders == null) {
            ArrayList<String> listOfHeaders = new ArrayList<String>();
            BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
            String line = null;
            while ((line = reader.readLine()) != null)
                if (line.startsWith(">")) {
                    if (truncateAfterFirstSpace && line.indexOf(' ') > -1)
                        listOfHeaders.add(line.substring(1, line.indexOf(' ')));
                    else
                        listOfHeaders.add(line.substring(1));
                }
            reader.close();

            if (storeCopy)
                this.listOfHeaders = (ArrayList<String>)listOfHeaders.clone();

            return listOfHeaders;
        }
        return (ArrayList<String>) this.listOfHeaders.clone();
    }

    @SuppressWarnings("unchecked")
    public HashMap<String, String> getHeaderToSequenceMap(boolean storeCopy) throws IOException {


        if (mapOfHeadersToSequences == null) {

            HashMap<String, String> mapOfHeadersToSequences = new HashMap<String, String>();
            BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
            String line;
            String headerCurrent = "";
            StringBuffer sequenceCurrent = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (headerCurrent.length() > 0)
                        mapOfHeadersToSequences.put(headerCurrent, sequenceCurrent.toString());
                    headerCurrent = line.substring(1).replaceAll("[\r\n]", "");
                    sequenceCurrent = new StringBuffer();
                } else
                    sequenceCurrent.append(line.trim().replaceAll("[\r\n]", ""));
            }
            if (headerCurrent.length() > 0)
                mapOfHeadersToSequences.put(headerCurrent, sequenceCurrent.toString());

            reader.close();

            if (storeCopy)
                this.mapOfHeadersToSequences = (HashMap<String, String>) mapOfHeadersToSequences.clone();

            return mapOfHeadersToSequences;
        }


        return (HashMap<String, String>) mapOfHeadersToSequences.clone();
    }

    public StringBuffer getSequenceByHeaderPrefix(String headerPrefixWithoutGTSign) throws Exception {

        return getSequenceByHeader(headerPrefixWithoutGTSign, true);
    }

    public StringBuffer getSequenceByHeader(String headerWithoutGTSign) throws Exception {

        return getSequenceByHeader(headerWithoutGTSign, false);
    }


    private StringBuffer getSequenceByHeader(String headerWithoutGTSign, boolean justRequirePrefix) throws IOException {

        StringBuffer sequence = null;

        BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
        String line = null;
        boolean sequenceHeaderFound = false;
        while (!sequenceHeaderFound && (line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if ((justRequirePrefix && line.startsWith(">" + headerWithoutGTSign)) ||
                        line.equalsIgnoreCase(">" + headerWithoutGTSign)) {
                    sequenceHeaderFound = true;
                    sequence = new StringBuffer();
                    while ((line = reader.readLine()) != null  && !line.startsWith(">"))
                        sequence.append(line.trim());

                }
            }
        }

        reader.close();

        return sequence;

    }


    public HashSet<String> writeSequencesByHeaders(File fileOutputFileFasta, HashSet<String> setOfHeadersWithoutGTSign, boolean keepNewlineFormatting, boolean truncateAfterFirstSpace) throws IOException {

        HashSet<String> setOfHeadersFound = new HashSet<String>();

        BufferedWriter writerOutputFileFasta = new BufferedWriter(new FileWriter(fileOutputFileFasta));
        BufferedReader readerInputFastaFile = new BufferedReader(new FileReader(this.fileFastaDatabase));
        String line = readerInputFastaFile.readLine();
        while (line != null) {
            String headerPossible = line.substring(1).trim();
            if (truncateAfterFirstSpace && headerPossible.indexOf(' ') > -1)
                headerPossible = headerPossible.substring(0 , headerPossible.indexOf(' '));
            if (setOfHeadersWithoutGTSign.contains(headerPossible)) {
                setOfHeadersFound.add(headerPossible);
                writerOutputFileFasta.write(line.trim());
                writerOutputFileFasta.newLine();
                while ((line = readerInputFastaFile.readLine()) != null  && !line.startsWith(">")) {
                    writerOutputFileFasta.write(line.trim());
                    if (keepNewlineFormatting)
                        writerOutputFileFasta.newLine();
                }
                if (!keepNewlineFormatting) writerOutputFileFasta.newLine();
            } else
                line = readerInputFastaFile.readLine();
        }

        readerInputFastaFile.close();

        HashSet<String> setOfHeadersNotFound = new HashSet<String>();
        Iterator<String> iteratorOverHeaders = setOfHeadersWithoutGTSign.iterator();
        while (iteratorOverHeaders.hasNext()) {
            String header = (String) iteratorOverHeaders.next();
            if (!setOfHeadersFound.contains(header))
                setOfHeadersNotFound.add(header);
        }

        return setOfHeadersNotFound;
    }

    public void generatePaddedVersionToFile(int lengthOfPadOnFivePrimeEnd, int lengthOfPadOnThreePrimeEnd, char characterToPadWith, File fileOutputFastaPadded) throws IOException {

        String sequenceFivePrimePad = "";
        for (int i = 0; i < lengthOfPadOnFivePrimeEnd; i++)
            sequenceFivePrimePad += characterToPadWith;
        String sequenceThreePrimePad = "";
        for (int i = 0; i < lengthOfPadOnThreePrimeEnd; i++)
            sequenceThreePrimePad += characterToPadWith;

        BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
        BufferedWriter writer = new BufferedWriter(new FileWriter(fileOutputFastaPadded));
        String line = null;
        String sequenceCurrent = "";
        String headerCurrent = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (headerCurrent.length() > 0) {
                    writer.write(headerCurrent);
                    writer.newLine();
                    writer.write(sequenceFivePrimePad + sequenceCurrent.toUpperCase() + sequenceThreePrimePad);
                    writer.newLine();
                }

                headerCurrent = line;
                sequenceCurrent = "";
            } else
                sequenceCurrent += line.trim();
        }
        if (headerCurrent.length() > 0) {
            writer.write(headerCurrent);
            writer.newLine();
            writer.write(sequenceFivePrimePad + sequenceCurrent.toUpperCase() + sequenceThreePrimePad);
            writer.newLine();
        }
        reader.close();
        writer.close();
    }

    /**
     * Truncate each sequence in a fasta file to a new fasta file.  Sequences are trimmed from the three prime end.
     *
     * @param lengthOfTruncatedReads
     * @param fileOutputFasta
     * @throws IOException
     */
    public void generateTruncatedVersionToFile(int lengthOfTruncatedReads,
                                                File fileOutputFasta) throws IOException {
        generateSplitVersionsToFiles(lengthOfTruncatedReads, -1, null, fileOutputFasta, null);

    }

    public void generateSplitVersionsToFiles(int lengthOfFirstPart, int lengthOfSecondPart,
                                             String stringToAppendToStartOfSecondPart,
                                             File fileOutputFastaFirstPart, File fileOutputFastaSecondPart) throws IOException {

        BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
        BufferedWriter writerFastaFirstPart = new BufferedWriter(new FileWriter(fileOutputFastaFirstPart));
        BufferedWriter writerFastaSecondPart = null;
        if (fileOutputFastaSecondPart != null)
            writerFastaSecondPart = new BufferedWriter(new FileWriter(fileOutputFastaSecondPart));
        String line = null;
        String sequenceCurrent = "";
        String headerCurrent = "";
        this.numberOfSequencesInDB = 0;
        Properties p = new Properties();
        p.setProperty(FileHeader.CREATED_BY_KEY, this.getClass().getSimpleName()+".generateSplitVersionsToFiles()");
        p.setProperty(FileHeader.SEQUENCE_SOURCE_FILE, this.fileFastaDatabase.getAbsolutePath());
        p.setProperty(FileHeader.SPLIT_KEY, "left");
        p.setProperty(FileHeader.SPLIT_LENGTH_KEY, Integer.toString(lengthOfFirstPart-1));
        FileHeader.writeHeader(p, writerFastaFirstPart, null);
        if (writerFastaSecondPart != null) {
        	p.setProperty(FileHeader.SPLIT_KEY, "right");
            p.setProperty(FileHeader.SPLIT_LENGTH_KEY, Integer.toString(lengthOfSecondPart));
            FileHeader.writeHeader(p, writerFastaSecondPart, null);
        }
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                this.numberOfSequencesInDB++;
                printSplitRead(lengthOfFirstPart, lengthOfSecondPart,
						stringToAppendToStartOfSecondPart,
						writerFastaFirstPart, writerFastaSecondPart,
						sequenceCurrent, headerCurrent);

                headerCurrent = line;
                sequenceCurrent = "";
            } else
                sequenceCurrent += line.trim();
        }
        if (headerCurrent.length() > 0)
            printSplitRead(lengthOfFirstPart, lengthOfSecondPart,
					stringToAppendToStartOfSecondPart,
					writerFastaFirstPart, writerFastaSecondPart,
					sequenceCurrent, headerCurrent);

        reader.close();
        writerFastaFirstPart.close();
        if (writerFastaSecondPart != null)
            writerFastaSecondPart.close();
    }

	private void printSplitRead(int lengthOfFirstPart, int lengthOfSecondPart,
			String stringToAppendToStartOfSecondPart,
			BufferedWriter writerFastaFirstPart,
			BufferedWriter writerFastaSecondPart, String sequenceCurrent,
			String headerCurrent) throws IOException {
		if (headerCurrent.length() > 0) {

		    writerFastaFirstPart.write(headerCurrent);
		    writerFastaFirstPart.newLine();
		    writerFastaFirstPart.write(sequenceCurrent.substring(0, lengthOfFirstPart).toUpperCase());
		    writerFastaFirstPart.newLine();

		    if (writerFastaSecondPart != null) {
		        writerFastaSecondPart.write(headerCurrent);
		        writerFastaSecondPart.newLine();
		        if (stringToAppendToStartOfSecondPart != null)
		            writerFastaSecondPart.write(stringToAppendToStartOfSecondPart);
		        writerFastaSecondPart.write(sequenceCurrent.substring(sequenceCurrent.length() - lengthOfSecondPart).toUpperCase());
		        writerFastaSecondPart.newLine();
		    }
		}
	}

    /**
     * FIXME: used a greedy approach.  splitting decisions made here could probably be optimized somehow.
     *
     * @param maxRefSequenceLength
     * @return
     * @throws IOException
     */
    public ArrayList<HashSet<String>> getSetsOfSequenceHeadersEachTotalingLessThanNMegabases(long maxRefSequenceLength) throws Exception {

        SortedMap<String, Long> mapHeaderToSequenceLength = getSequenceLengths();

        this.numberOfSequencesInDB = mapHeaderToSequenceLength.keySet().size();

        ArrayList<HashSet<String>> listOfHeaderSets = new ArrayList<HashSet<String>>();
        ArrayList<Long> listOfHeaderSetSequenceLengths = new ArrayList<Long>();

        Iterator<String> iteratorOverHeaders = mapHeaderToSequenceLength.keySet().iterator();
        while (iteratorOverHeaders.hasNext()) {
            String headerOfCurrentSequence = (String) iteratorOverHeaders.next();
            long lengthOfCurrentSequence = mapHeaderToSequenceLength.get(headerOfCurrentSequence);

            if (lengthOfCurrentSequence > maxRefSequenceLength)
                logger.info("WARNING: Reference sequence " + headerOfCurrentSequence + " has length " + lengthOfCurrentSequence + ", which is longer than the max reference sequence length specified (" + maxRefSequenceLength + " ).");

            boolean foundExisitingHeaderSetWithRoom = false;
            int indexOfCurrentHeaderSet = -1;
            while (!foundExisitingHeaderSetWithRoom && ++indexOfCurrentHeaderSet < listOfHeaderSets.size()) {
                if (listOfHeaderSetSequenceLengths.get(indexOfCurrentHeaderSet) + lengthOfCurrentSequence <=  maxRefSequenceLength) {
                    foundExisitingHeaderSetWithRoom = true;
                    listOfHeaderSets.get(indexOfCurrentHeaderSet).add(headerOfCurrentSequence);
                    listOfHeaderSetSequenceLengths.set(indexOfCurrentHeaderSet, listOfHeaderSetSequenceLengths.get(indexOfCurrentHeaderSet) + lengthOfCurrentSequence);
                }
            }
            if (!foundExisitingHeaderSetWithRoom) {
                HashSet<String> setOfHeadersNew = new HashSet<String>();
                setOfHeadersNew.add(headerOfCurrentSequence);
                listOfHeaderSets.add(setOfHeadersNew);
                listOfHeaderSetSequenceLengths.add(lengthOfCurrentSequence);
            }

        }        

        return listOfHeaderSets;
    }


    public SortedMap<String, Long> getSequenceLengths() throws Exception {
    	return getSequenceLengths(-1);
    }
    
    /**
     * 
     * @return Map of Sequence Identifiers to sequence lengths sorted by identifiers
     * @throws IOException
     */
    public SortedMap<String, Long> getSequenceLengths(long max) throws Exception {
        TreeMap<String, Long> mapHeaderToSequenceLength = new TreeMap<String, Long>();

        ArrayList<String> listOfHeaders = new ArrayList<String>();

        BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
        String line = null;
        long lengthOfCurrentSequence = 0;
        String headerOfCurrentSequence = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
            	if (max > -1 && mapHeaderToSequenceLength.size() > max) break;
                if (headerOfCurrentSequence.length() > 0)
                    mapHeaderToSequenceLength.put(headerOfCurrentSequence, lengthOfCurrentSequence);
                headerOfCurrentSequence = line.substring(1).trim();
                lengthOfCurrentSequence = 0;
                listOfHeaders.add(headerOfCurrentSequence);
            } else
                lengthOfCurrentSequence += line.trim().length();
        }
        if (headerOfCurrentSequence.length() > 0)
            mapHeaderToSequenceLength.put(headerOfCurrentSequence, lengthOfCurrentSequence);

        reader.close();

        if (max < 0)  { //If we looked at every sequence, update the instance variables.
        	this.listOfHeaders = listOfHeaders;
        	this.numberOfSequencesInDB = mapHeaderToSequenceLength.keySet().size();
        }

        return mapHeaderToSequenceLength;
    }
    
    public Long getMonoLength(int numberOfSeqsToCheck) throws Exception {
    	Long monoLength = null;
    	for (Long value : getSequenceLengths(numberOfSeqsToCheck).values()) {
    		if (monoLength == null) {
    			monoLength = value;
    			continue;
    		}
    		if (!value.equals(monoLength)) throw new Exception("This file contains sequences of varying length.");
    	}
    	return monoLength;
    }


    public long getTotalLengthOfSequence() throws Exception {

        long lengthOfAllSequence = 0;

        BufferedReader reader = new BufferedReader(new FileReader(this.fileFastaDatabase));
        String line;
        while ((line = reader.readLine()) != null)
            if (!line.startsWith(">"))
                lengthOfAllSequence += line.trim().length();
        
        reader.close();

        return lengthOfAllSequence;

    }
    
    public static List<Integer>[] partitionFastaFileToFastaFilesByHeaderSets(File fileReferencePartitionTable,
            File[] filesReferenceFastaPartitions,
            FastaDatabase fastaDBOfFullReference,
            ArrayList<HashSet<String>> listOfReferenceHeaderSetsForEachPartition) throws Exception {
    	PrintWriter writerReferencePartitionTable = null;
    	try {
    		writerReferencePartitionTable = new PrintWriter(new BufferedWriter(new FileWriter(fileReferencePartitionTable)));
    		return partitionFastaFileToFastaFilesByHeaderSets(writerReferencePartitionTable, filesReferenceFastaPartitions, fastaDBOfFullReference, listOfReferenceHeaderSetsForEachPartition);
    	} finally {
    		if (writerReferencePartitionTable != null) writerReferencePartitionTable.close();
    	}
    }
    
    public static List<Integer>[] partitionFastaFileToFastaFilesByHeaderSets(PrintWriter writerReferencePartitionTable,
                                                                             File[] filesReferenceFastaPartitions,
                                                                             FastaDatabase fastaDBOfFullReference,
                                                                             List<HashSet<String>> listOfReferenceHeaderSetsForEachPartition) throws Exception {
        //ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
        List<String> listOfReferenceHeadersOrdered = fastaDBOfFullReference.getListOfHeaders(false, false);
        HashMap<String, Integer> mapHeaderToReferenceSequenceNumber = new HashMap<String, Integer>();
        for (int i = 0; i < listOfReferenceHeadersOrdered.size(); i++)
             mapHeaderToReferenceSequenceNumber.put(listOfReferenceHeadersOrdered.get(i), i +1);
        
        writerReferencePartitionTable.println("#refPartitionIndex\tsequenceNumberInRefPartition(1-based)\tsequenceNumberInFullReference(1-based)\tsequenceHeaderInFullReference");
        
        @SuppressWarnings("unchecked")
        ArrayList<Integer>[] listsOfFullRefSequenceNumbersOrderedForEachReferencePartition = new ArrayList[listOfReferenceHeaderSetsForEachPartition.size()];
        
        for (int indexPartition = 0; indexPartition < filesReferenceFastaPartitions.length; indexPartition++) {
            //filesReferenceFastaPartitions[indexPartition] = new File(folderForTempFiles.getPath() + "/reference." + indexPartition + ".fasta");
            ArrayList<Integer> listOfRefSequencesOrdered = new ArrayList<Integer>();
            logger.info(" -->Creating split reference file " + indexPartition + "...");
            BufferedWriter writerFastaFilePartition = new BufferedWriter(new FileWriter(filesReferenceFastaPartitions[indexPartition]));
            /*Iterator iteratorOverHeaderSet = listOfReferenceHeaderSetsForEachPartition.get(indexPartition).iterator();
            while (iteratorOverHeaderSet.hasNext()) {
                String header = (String) iteratorOverHeaderSet.next();
                writerFastaFilePartition.write(">" + header);
                writerFastaFilePartition.newLine();
                // FIXME: this is very slow when the number of reference sequences is high
                fastaDBOfFullReference.writeSequenceByHeaderPrefix(writerFastaFilePartition, header, false);
                writerFastaFilePartition.newLine();

                int numberOfSequenceInFullReference = mapHeaderToReferenceSequenceNumber.get(header);
                listOfRefSequencesOrdered.add(numberOfSequenceInFullReference);
                writerReferencePartitionTable.write(indexPartition + "\t" + listOfRefSequencesOrdered.size() + "\t" + numberOfSequenceInFullReference + "\t" + header);
                writerReferencePartitionTable.newLine();

            } */

            HashSet<String> headersForThisPartition = listOfReferenceHeaderSetsForEachPartition.get(indexPartition);
            BufferedReader readerInputFastaFile = new BufferedReader(new FileReader(fastaDBOfFullReference.fileFastaDatabase));
            String line = readerInputFastaFile.readLine();
            while (line != null) {
                if (line.startsWith(">") && headersForThisPartition.contains(line.substring(1).trim())) {
                    String header = line.substring(1).trim();
                    writerFastaFilePartition.write(">" + header);
                    writerFastaFilePartition.newLine();
                    while ((line = readerInputFastaFile.readLine()) != null  && !line.startsWith(">")) {
                        writerFastaFilePartition.write(line.trim());
                        //writerOutputFileFasta.newLine();
                    }
                    writerFastaFilePartition.newLine();

                    int numberOfSequenceInFullReference = mapHeaderToReferenceSequenceNumber.get(header);
                    listOfRefSequencesOrdered.add(numberOfSequenceInFullReference);
                    writerReferencePartitionTable.println(indexPartition + "\t" + listOfRefSequencesOrdered.size() + "\t" + numberOfSequenceInFullReference + "\t" + header);

                } else
                    line = readerInputFastaFile.readLine();
            }


            writerFastaFilePartition.close();

            //fastaDBOfFullReference.writeSequencesByHeaders(filesReferenceFastaPartitions[indexPartition], listOfReferenceHeaderSetsForEachPartition.get(indexPartition), false);

            listsOfFullRefSequenceNumbersOrderedForEachReferencePartition[indexPartition] = listOfRefSequencesOrdered;
        }
        writerReferencePartitionTable.close();
        return listsOfFullRefSequenceNumbersOrderedForEachReferencePartition;
    }

	public File getFileFastaDatabase() {
		return fileFastaDatabase;
	}


}
