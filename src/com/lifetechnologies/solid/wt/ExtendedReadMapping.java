package com.lifetechnologies.solid.wt;

import java.util.*;
import java.io.*;

/**
 * User: tuchbb
 * Date: Oct 1, 2008
 * Time: 11:28:17 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExtendedReadMapping extends ReadMapping {

    private int score;
    private int positionOfAlignmentStartInRead;
    private int lengthOfAlignment;
    private String sequenceOfRead;
    private TreeSet<Integer> setOfStartPositionsForValidAdjacentMismatches = new TreeSet<Integer>();
    private ArrayList<Insertion> listOfInsertionsInRead = new ArrayList<Insertion>();
    private ArrayList<Insertion> listOfInsertionsInReference = new ArrayList<Insertion>();

    public ExtendedReadMapping(String idOfRead, String idOfMatchingReferenceSequence, int indexOfMatchingReferenceSequence, int positionOfAlignmentStartInReferenceSequence, int numberOfMismatches, int score, int positionOfAlignmentStartInRead, int lengthOfAlignment) {
        super(idOfRead, idOfMatchingReferenceSequence, indexOfMatchingReferenceSequence, positionOfAlignmentStartInReferenceSequence, numberOfMismatches);
        this.score = score;
        this.positionOfAlignmentStartInRead = positionOfAlignmentStartInRead;
        this.lengthOfAlignment = lengthOfAlignment;
    }

    public String getSequenceOfRead() {
        return sequenceOfRead;
    }

    public void setSequenceOfRead(String sequenceOfRead) {
        this.sequenceOfRead = sequenceOfRead;
    }

    public void addValidAdjacenStartPosition(int position) {
        setOfStartPositionsForValidAdjacentMismatches.add(position);
    }

    public void addInsertionInRead(Insertion insertion) {
        listOfInsertionsInRead.add(insertion);
    }

    public void addInsertionInReference(Insertion insertion) {
        listOfInsertionsInReference.add(insertion);
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public int getPositionOfAlignmentStartInRead() {
        return positionOfAlignmentStartInRead;
    }

    public void setPositionOfAlignmentStartInRead(int positionOfAlignmentStartInRead) {
        this.positionOfAlignmentStartInRead = positionOfAlignmentStartInRead;
    }

    public int getPositionOfAlignmentEndInRead() {
        return positionOfAlignmentStartInRead + this.lengthOfAlignment -1;
    }

    public int getLengthOfAlignment() {
        return lengthOfAlignment;
    }

    public void setLengthOfAlignment(int lengthOfAlignment) {
        this.lengthOfAlignment = lengthOfAlignment;
    }

    public TreeSet<Integer> getSetOfStartPositionsForValidAdjacentMismatches() {
        return setOfStartPositionsForValidAdjacentMismatches;
    }

    public void setSetOfStartPositionsForValidAdjacentMismatches(TreeSet<Integer> setOfStartPositionsForValidAdjacentMismatches) {
        this.setOfStartPositionsForValidAdjacentMismatches = setOfStartPositionsForValidAdjacentMismatches;
    }

    public ArrayList<Insertion> getListOfInsertionsInRead() {
        return listOfInsertionsInRead;
    }

    public void setListOfInsertionsInRead(ArrayList<Insertion> listOfInsertionsInRead) {
        this.listOfInsertionsInRead = listOfInsertionsInRead;
    }

    public ArrayList<Insertion> getListOfInsertionsInReference() {
        return listOfInsertionsInReference;
    }

    public void setListOfInsertionsInReference(ArrayList<Insertion> listOfInsertionsInReference) {
        this.listOfInsertionsInReference = listOfInsertionsInReference;
    }

    /**
     * This is the position of the base in the reference that aligns with
     * the first color, regardless of whether it matches.
     */ 
    public GenomicPosition getPositionOfFirstColor() {
    	if (getPositionOfAlignmentStartInReferenceSequence() > 0)
    		return new GenomicPosition(
    				this.getIndexOfMatchingReferenceSequence(),
    				getPositionOfAlignmentStartInReferenceSequence() - this.getPositionOfAlignmentStartInRead() + 1,
    				Strand.POSITIVE);
    	else
    		return new GenomicPosition(
    				this.getIndexOfMatchingReferenceSequence(),
    				Math.abs(getPositionOfAlignmentStartInReferenceSequence()) + this.getLengthOfAlignment() + this.getPositionOfAlignmentStartInRead() - 2,
    				Strand.NEGATIVE);
    }
    
    /**
     * This is the postition of the base in the reference that aligns with
     * the last color, regardless of whether it matches.
     * @param readLength
     * @return
     */
    public GenomicPosition getPositionOfLastColor(int readLength) {
    	if (readLength < 1) throw new IllegalArgumentException("readLength cannot be less than 1");
    	if (getPositionOfAlignmentStartInReferenceSequence() > 0)
    		return getPositionOfFirstColor().derivePosition(readLength - 1);
    	else
    		return getPositionOfFirstColor().derivePosition(1 - readLength);
    }

	public String toString() {

        StringBuffer header = new StringBuffer();

        // append basic info
        header.append(",");
        header.append(this.getIdOfMatchingReferenceSequence());
        header.append(".");
        header.append(this.getPositionOfAlignmentStartInReferenceSequence());
        header.append(".");
        header.append(this.getPositionOfAlignmentStartInRead());
        header.append(".");
        header.append(this.getLengthOfAlignment());
        header.append(".");
        header.append(this.getScore());
        header.append(".");
        header.append(this.getNumberOfMismatches());

        // append valid adjacent positions
        header.append("#");

        //Valid adjacent reporting Disabled by default.
        if (Utilities.isTrue(System.getProperty(WholeTranscriptomeAnalyzer.EXPERIMENTAL_MODE_SYSTEM_PROPERTY))) {
        	Iterator<Integer> iteratorOverValidAdjacents = this.getSetOfStartPositionsForValidAdjacentMismatches().iterator();
        	while (iteratorOverValidAdjacents.hasNext()) {
        		header.append(iteratorOverValidAdjacents.next());
        		header.append(".");
        	}
 
        	if (this.getSetOfStartPositionsForValidAdjacentMismatches().size() > 0)
        		header.deleteCharAt(header.length() -1);
        }
        // append insertions in read
        header.append("#");
        Iterator<Insertion> iteratorOverInsertions = this.getListOfInsertionsInRead().iterator();
        while (iteratorOverInsertions.hasNext()) {
            header.append(iteratorOverInsertions.next() );
            header.append(";");
        }
        if (this.getListOfInsertionsInRead().size() > 0)
            header.deleteCharAt(header.length() -1);

        // append insertions in reference
        header.append("#");
        iteratorOverInsertions = this.getListOfInsertionsInReference().iterator();
        while (iteratorOverInsertions.hasNext()) {
            header.append(iteratorOverInsertions.next() );
            header.append(";");
        }
        if (this.getListOfInsertionsInReference().size() > 0)
            header.deleteCharAt(header.length() -1);


        return header.toString();
    }


    public static void main(String args[]) throws IOException {
        String line = null;
        int count =0;
        long charsTOSKip = 3469675108L;
        long startTime = System.currentTimeMillis();
        int numberOfPieces = (int)Math.ceil((double)charsTOSKip / 2000000000);

        BufferedReader readerA = new BufferedReader(new FileReader("/home/tuchbb/projects/MAQC/mapping_results/73_20080930_1_Ambion_HBR_HBR_F3.csfasta"));
        startTime = System.currentTimeMillis();
        while ((line = readerA.readLine()) != null) {
            count++;
            if (count == 50000000) {
                System.out.println(line);
                System.out.println("Time A1\t" + (System.currentTimeMillis() - startTime) / 1000.0);
                break;
            }
            readerA.readLine();
        }

        startTime = System.currentTimeMillis();
        count =0;
        while (count++ < 5000000)
            line = readerA.readLine();
        readerA.close();
        System.out.println(line);
        System.out.println("Time A2\t" + (System.currentTimeMillis() - startTime) / 1000.0);


        BufferedReader readerB = new BufferedReader(new FileReader("/home/tuchbb/projects/MAQC/mapping_results/73_20080930_1_Ambion_HBR_HBR_F3.csfasta"));
        for (int i = 0; i < numberOfPieces -1; i++) 
            readerB.skip(2000000000);
        readerB.skip(charsTOSKip % 2000000000);
        System.out.println(readerB.readLine());
        System.out.println("Time B1\t" + (System.currentTimeMillis() - startTime) / 1000.0);

        startTime = System.currentTimeMillis();
        count =0;
        while (count++ < 5000000)
            line = readerB.readLine();
        readerB.close();
        System.out.println(line);
        System.out.println("Time B2\t" + (System.currentTimeMillis() - startTime) / 1000.0);


        RandomAccessFile randomAccessFile = new RandomAccessFile("/home/tuchbb/projects/MAQC/mapping_results/73_20080930_1_Ambion_HBR_HBR_F3.csfasta", "r");
        startTime = System.currentTimeMillis();
        for (int i = 0; i < numberOfPieces -1; i++)
            randomAccessFile.skipBytes(2000000000);
        randomAccessFile.skipBytes((int)(charsTOSKip % 2000000000));
        System.out.println(randomAccessFile.readLine());
        System.out.println("Time C1\t" + (System.currentTimeMillis() - startTime) / 1000.0);

        startTime = System.currentTimeMillis();
        count =0;
        while (count++ < 5000000)
            line = randomAccessFile.readLine();
        System.out.println(line);
        System.out.println("Time C2\t" + (System.currentTimeMillis() - startTime) / 1000.0);
        randomAccessFile.close();



        startTime = System.currentTimeMillis();
        BufferedRandomAccessFile bufferedRandomAccessFile = new BufferedRandomAccessFile(new File("/home/tuchbb/projects/MAQC/mapping_results/73_20080930_1_Ambion_HBR_HBR_F3.csfasta"), "r", 1024 * 64);
        for (int i = 0; i < numberOfPieces -1; i++)
            bufferedRandomAccessFile.skipBytes(2000000000);
        bufferedRandomAccessFile.skipBytes((int)(charsTOSKip % 2000000000));
        System.out.println(bufferedRandomAccessFile.readLine());
        System.out.println("Time D1\t" + (System.currentTimeMillis() - startTime) / 1000.0);

        startTime = System.currentTimeMillis();
        count =0;
        while (count++ < 5000000)
            line = bufferedRandomAccessFile.readNextLine();
        System.out.println(line);
        System.out.println("Time D2\t" + (System.currentTimeMillis() - startTime) / 1000.0);
        bufferedRandomAccessFile.close();

    }

}
