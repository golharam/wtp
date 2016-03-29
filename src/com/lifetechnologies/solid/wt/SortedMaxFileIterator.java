package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

import com.lifetechnologies.util.SeekableIterator;

/**
 * Seek-able and Peek-able Iteration through a Max File.
 * @author mullermw
 *
 */
public class SortedMaxFileIterator implements SeekableIterator<ExtendedReadMapping> {

	private BufferedRandomAccessFile sortedMaxFile;
	private ExtendedReadMappingEntry firstEntry, lastEntry, nextEntry;
	private long lengthOfReads = -1;
	
	public static final Comparator<ExtendedReadMapping> REFERENCE_POSITION_COMPARATOR = new Comparator<ExtendedReadMapping>() {
		@Override
		public int compare(ExtendedReadMapping o1, ExtendedReadMapping o2) {
			if (o1 == null ) {
				if (o2 == null) return 0;
				return -1;
			}
			if (o2 == null) return 1;
			if (o1.equals(o2)) return 0;
			if (o1.getIndexOfMatchingReferenceSequence() != o2.getIndexOfMatchingReferenceSequence())
				return o1.getIndexOfMatchingReferenceSequence() - o2.getIndexOfMatchingReferenceSequence();
			if (o1.getPositionOfAlignmentStartInReferenceSequence() != o2.getPositionOfAlignmentStartInReferenceSequence())
				return Math.abs(o1.getPositionOfAlignmentStartInReferenceSequence()) - Math.abs(o2.getPositionOfAlignmentStartInReferenceSequence());
			return 0;
		}	
	};

	public SortedMaxFileIterator(File _sortedMaxFile) throws IOException {
		sortedMaxFile = new BufferedRandomAccessFile(_sortedMaxFile, "r");
		for (String line=sortedMaxFile.readNextLine(); line != null; line = sortedMaxFile.readNextLine() ) {
			if (line.startsWith(">") == false) continue;
			long offset = sortedMaxFile.getFilePointer() - line.length() - 1;
			line = line.substring(1);
			firstEntry = new ExtendedReadMappingEntry(line, offset);
			line = sortedMaxFile.readNextLine();
			lengthOfReads = line.length();
			break;
		}

		sortedMaxFile.seek(sortedMaxFile.length() - 1000);

		for (String line=sortedMaxFile.readNextLine(); line != null; line = sortedMaxFile.readNextLine() ) {
			if (line.startsWith(">") == false) continue;
			long offset = sortedMaxFile.getFilePointer() - line.length() - 1;
			line = line.substring(1);
			lastEntry = new ExtendedReadMappingEntry(line, offset);
		}
		sortedMaxFile.seek(0);
	}

	@Override
	public boolean hasNext() {
		try {
			return readNextEntry(true) != null;
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	@Override
	public ExtendedReadMapping next() {
		try {
			ExtendedReadMapping e = readNextEntry(true);
			readNextEntry(false);
			return e;
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public void remove() {
		throw new UnsupportedOperationException("remove() is not supported");
	}

	public void seekTo(ExtendedReadMapping prototype) {
		try {
			//Seeking to before the first mapping.
			if (REFERENCE_POSITION_COMPARATOR.compare(prototype, firstEntry) <= 0) {
				sortedMaxFile.seek(firstEntry.getOffset());
				return;
			}
			
			//Seeking beyond the last mapping.
			if (REFERENCE_POSITION_COMPARATOR.compare(prototype, lastEntry) > 0) {
				sortedMaxFile.seek(sortedMaxFile.length()-1);
				return;
			}
			
			//binary search for a read(entry) that is greater than or equal to 
			//prototype and preceded by one that is less than the prototype
			long position = -1;
			long left = 0;
			long right = sortedMaxFile.length() - 1;
			while (position < 0) {
				long i = (left + right) / 2;
				sortedMaxFile.seek(i);
				ExtendedReadMappingEntry precedingEntry = readNextEntry(false);
				ExtendedReadMappingEntry entry = readNextEntry(false);
				//System.out.println(precedingEntry + " " + entry + " " + REFERENCE_POSITION_COMPARATOR.compare(prototype, precedingEntry) + " " + REFERENCE_POSITION_COMPARATOR.compare(prototype, entry));
				if (REFERENCE_POSITION_COMPARATOR.compare(prototype, precedingEntry) <= 0) {
					right = i;
				} else {
					if (REFERENCE_POSITION_COMPARATOR.compare(prototype, entry) <= 0) {
						//Got it, break out of loop.
						position = entry.getOffset();
					} else {
						left = i;
					}
				}
			}
			return;
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	public ExtendedReadMapping peek() {
		try {
			if (readNextEntry(true) != null)
				return nextEntry.mapping;
			return null;
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	

	public long getOffsetOfNextEntry() {
		if (nextEntry == null) return -1;
		return nextEntry.getOffset();
	}
	
	public static ExtendedReadMapping createPrototype(int refidx, int pos) {
		return new ExtendedReadMapping("", "", refidx, pos, 0, 0, 0, 0);
	}

	private ExtendedReadMappingEntry readNextEntry(boolean readFromBuffer) throws IOException {
		if (nextEntry == null || !readFromBuffer) {
			nextEntry = null;
			String line;
			for (line=sortedMaxFile.readNextLine(); line != null; line = sortedMaxFile.readNextLine() ) {
				//System.out.println(line);
				if (line.startsWith(">") == false) {
					if (nextEntry == null)
						continue;
					else {
						nextEntry.setSequenceOfRead(line);
						break;
					}
				}
				long offset = sortedMaxFile.getFilePointer() - line.length() - 1;
				nextEntry = new ExtendedReadMappingEntry(line, offset);
			}
			if (line == null) nextEntry = null;
		}
		return nextEntry;
	}

	@Override
	protected void finalize() throws Throwable {
		sortedMaxFile.close();
		super.finalize();
	}

	public static void main(String[] args) {
		try {
			SortedMaxFileIterator iter = new SortedMaxFileIterator(new File(args[0]));

			RandomAccessFile file = new RandomAccessFile(new File(args[0]), "r");
			ExtendedReadMapping prototype = createPrototype(Integer.parseInt(args[1]), Integer.parseInt(args[2]));
			iter.seekTo(prototype);
			file.seek(iter.getOffsetOfNextEntry());
			for (int i=0; i<10; i++) {
				System.out.println(iter.next());
			}
		} catch (Exception e ) {
			e.printStackTrace();
		}
	}

	public long getLengthOfReads() {
		return lengthOfReads;
	}

	public void setLengthOfReads(long lengthOfReads) {
		this.lengthOfReads = lengthOfReads;
	}
}
	
class ExtendedReadMappingEntry extends ExtendedReadMapping {
	long offset = -1;
	ExtendedReadMapping mapping;

	public long getOffset() {
		return offset;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

	public ExtendedReadMappingEntry(String line, long offset) {
		super(null, null, -1, -1, -1, -1, -1, -1);
		mapping = ExtendedReadMappingSetForRead.parseExtendedReadMappingsFromHeader(line)[0];
		this.offset = offset;
	}

	public void addInsertionInRead(Insertion insertion) {
		mapping.addInsertionInRead(insertion);
	}

	public void addInsertionInReference(Insertion insertion) {
		mapping.addInsertionInReference(insertion);
	}

	public void addValidAdjacenStartPosition(int position) {
		mapping.addValidAdjacenStartPosition(position);
	}

	public boolean equals(Object obj) {
		return mapping.equals(obj);
	}

	public String getIdOfMatchingReferenceSequence() {
		return mapping.getIdOfMatchingReferenceSequence();
	}

	public String getIdOfRead() {
		return mapping.getIdOfRead();
	}

	public int getIndexOfMatchingReferenceSequence() {
		return mapping.getIndexOfMatchingReferenceSequence();
	}

	public int getLengthOfAlignment() {
		return mapping.getLengthOfAlignment();
	}

	public ArrayList<Insertion> getListOfInsertionsInRead() {
		return mapping.getListOfInsertionsInRead();
	}

	public ArrayList<Insertion> getListOfInsertionsInReference() {
		return mapping.getListOfInsertionsInReference();
	}

	public int getNumberOfMismatches() {
		return mapping.getNumberOfMismatches();
	}

	public int getPositionOfAlignmentEndInRead() {
		return mapping.getPositionOfAlignmentEndInRead();
	}

	public int getPositionOfAlignmentStartInRead() {
		return mapping.getPositionOfAlignmentStartInRead();
	}

	public int getPositionOfAlignmentStartInReferenceSequence() {
		return mapping.getPositionOfAlignmentStartInReferenceSequence();
	}

	public int getScore() {
		return mapping.getScore();
	}

	public String getSequenceOfRead() {
		return mapping.getSequenceOfRead();
	}

	public TreeSet<Integer> getSetOfStartPositionsForValidAdjacentMismatches() {
		return mapping.getSetOfStartPositionsForValidAdjacentMismatches();
	}

	public int hashCode() {
		return mapping.hashCode();
	}

	public void setIdOfMatchingReferenceSequence(
			String idOfMatchingReferenceSequence) {
		mapping.setIdOfMatchingReferenceSequence(idOfMatchingReferenceSequence);
	}

	public void setIdOfRead(String idOfRead) {
		mapping.setIdOfRead(idOfRead);
	}

	public void setIndexOfMatchingReferenceSequence(
			int indexOfMatchingReferenceSequence) {
		mapping
				.setIndexOfMatchingReferenceSequence(indexOfMatchingReferenceSequence);
	}

	public void setLengthOfAlignment(int lengthOfAlignment) {
		mapping.setLengthOfAlignment(lengthOfAlignment);
	}

	public void setListOfInsertionsInRead(
			ArrayList<Insertion> listOfInsertionsInRead) {
		mapping.setListOfInsertionsInRead(listOfInsertionsInRead);
	}

	public void setListOfInsertionsInReference(
			ArrayList<Insertion> listOfInsertionsInReference) {
		mapping.setListOfInsertionsInReference(listOfInsertionsInReference);
	}

	public void setNumberOfMismatches(int numberOfMismatches) {
		mapping.setNumberOfMismatches(numberOfMismatches);
	}

	public void setPositionOfAlignmentStartInRead(
			int positionOfAlignmentStartInRead) {
		mapping
				.setPositionOfAlignmentStartInRead(positionOfAlignmentStartInRead);
	}

	public void setPositionOfAlignmentStartInReferenceSequence(
			int positionOfAlignmentStartInReferenceSequence) {
		mapping
				.setPositionOfAlignmentStartInReferenceSequence(positionOfAlignmentStartInReferenceSequence);
	}

	public void setScore(int score) {
		mapping.setScore(score);
	}

	public void setSequenceOfRead(String sequenceOfRead) {
		mapping.setSequenceOfRead(sequenceOfRead);
	}

	public void setSetOfStartPositionsForValidAdjacentMismatches(
			TreeSet<Integer> setOfStartPositionsForValidAdjacentMismatches) {
		mapping
				.setSetOfStartPositionsForValidAdjacentMismatches(setOfStartPositionsForValidAdjacentMismatches);
	}

	public String toString() {
		return mapping.toString();
	}

}
