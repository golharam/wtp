package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.lifetechnologies.solid.wt.FastaDatabaseIndex.SequenceInfo;
import com.lifetechnologies.util.BaseStream;

/**
 * Extends FastaDatabase to make use of an index.
 * @author mullermw
 *
 */
public class IndexedFastaDatabase extends FastaDatabase {

	private FastaDatabaseIndex index;
	private RandomAccessFile fastaFile;
	private File indexFile;
	private int indexOfNextSequence = 0;
	private boolean writeIndex = true;

	/**
	 * 
	 * 
	 * @param fastaFile the file containing the sequences.
	 * @param indexFile the index file
	 * @throws IOException if an IO problem occurs in reading the fasta or index file.
	 * @throws ParseException if there is a problem parsing the index or fasta file.
	 * @throws InvalidIndexException if the index is invalid.
	 */
	public IndexedFastaDatabase(File fastaFile, File indexFile, boolean writeIndex) throws IOException, ParseException, InvalidIndexException { 
		super();
		this.writeIndex = writeIndex;
		if (fastaFile == null) throw new IllegalArgumentException("fastaFile cannot be null.");
		if (indexFile == null) throw new IllegalArgumentException("indexFile cannot be null.");
		if (!fastaFile.exists()) throw new IllegalArgumentException(fastaFile + " doesn't exist.");
		this.indexFile = indexFile;
		this.fileFastaDatabase = fastaFile;
		if (indexFile.exists()) {
			FileInputStream fis = new FileInputStream(indexFile);
			try {
				this.index = FastaDatabaseIndex.newIndexFromIndexStream(fis);
				if (!this.index.getTarget().getAbsoluteFile().equals(fastaFile.getAbsoluteFile()))
						this.index = null;
			} finally {
				fis.close();
			}
		}
		if (this.index == null) {
			this.index = FastaDatabaseIndex.newIndexFromFasta(fastaFile);
			if (writeIndex) index.writeIndex(indexFile);
		}
		checkAndUpdateIndex();
		this.fastaFile = new RandomAccessFile(index.getTarget(), "r");
	}
	
	public IndexedFastaDatabase(File fastaFile, File indexFile) throws IOException, ParseException, InvalidIndexException { 
		this(fastaFile, indexFile, true);
	}
	
	private void checkAndUpdateIndex() throws IOException, ParseException, InvalidIndexException {
		if (index.isAccurate()) return;
		if (writeIndex) {
			index.rebuildIndex();
			index.writeIndex(indexFile);
		}
		if (!index.isAccurate()) throw new InvalidIndexException("Index is not accurate for " + index.getTarget(), this.index);
	}

		
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		if (fastaFile != null) fastaFile.close();
	}
	
	@Override
	public void reset() throws IOException {
		indexOfNextSequence = 0;
	}
	
	@Override
    public SequenceHeaderPair nextSequence(boolean toUpperCase) throws Exception {
		checkAndUpdateIndex();
		if (indexOfNextSequence > index.getNumberOfSequences() - 1) return null;
		FastaDatabaseIndex.SequenceInfo info = index.getSequenceInfo(indexOfNextSequence++);
		this.fastaFile.seek(info.getOffset());
        SequenceHeaderPair sequenceHeaderPair = null;
        StringBuffer stringBufferSequence = new StringBuffer();
        String line = this.fastaFile.readLine(); //Skip the header line.
        if (line == null) throw new Exception("The index file does not match the fasta file.  No data at file offset: "+info.getOffset());
        if (!line.startsWith(">")) throw new Exception("The index file does not match the fasta file.  Character at offset="+info.getOffset()+" is not a '>'");
        while ((line = this.fastaFile.readLine()) != null && !line.startsWith(">"))
            stringBufferSequence.append(line.trim());
        sequenceHeaderPair = new SequenceHeaderPair(">"+info.getHeader(), stringBufferSequence.toString());
        if (toUpperCase)    sequenceHeaderPair.setSequence(sequenceHeaderPair.getSequence().toUpperCase());
        return sequenceHeaderPair;
    }
	
	public BaseStream nextSequenceAsStream() throws IOException, InvalidIndexException, ParseException {
		checkAndUpdateIndex();
		if (indexOfNextSequence > index.getNumberOfSequences() - 1) return null;
		FastaDatabaseIndex.SequenceInfo info = index.getSequenceInfo(indexOfNextSequence++);
		final BufferedRandomAccessFile file = new BufferedRandomAccessFile(fileFastaDatabase, "r");
		file.seek(info.getOffset());
		file.readLine();
		return new MyBaseStream(file);
	}
	
	@Override
	public int getNumberOfSequencesInDatabase() throws IOException, ParseException, InvalidIndexException {
		checkAndUpdateIndex();
        return index.getNumberOfSequences();
    }
	
	@Override
	@SuppressWarnings("unchecked")
    public List<String> getListOfHeaders(boolean storeCopy, boolean truncateAfterFirstSpace) throws IOException, ParseException, InvalidIndexException {
		checkAndUpdateIndex();
		List<String> headers = new ArrayList<String>();
		for (FastaDatabaseIndex.SequenceInfo info : index.getSequenceInfo()) {
			if (truncateAfterFirstSpace) {
				headers.add(info.getSeqId());
			} else {
				headers.add(info.getHeader());
			}
		}
		return headers;
	}
	
	@Override
    public HashMap<String, String> getHeaderToSequenceMap(boolean storeCopy) throws IOException {
		return super.getHeaderToSequenceMap(storeCopy);
	}
	
	@Override
	public StringBuffer getSequenceByHeaderPrefix(
			String headerPrefixWithoutGTSign) throws IOException, InvalidIndexException, ParseException {
		return this.getSequenceByHeader(headerPrefixWithoutGTSign, true);
	}
	
	@Override
	public StringBuffer getSequenceByHeader(String headerWithoutGTSign)
			throws IOException, InvalidIndexException, ParseException {
		return this.getSequenceByHeader(headerWithoutGTSign, false);
	}
	
	private StringBuffer getSequenceByHeader(String header, boolean prefixOnly) throws IOException, InvalidIndexException, ParseException {
		checkAndUpdateIndex();
		FastaDatabaseIndex.SequenceInfo info = null;
		for (SequenceInfo _i : index.getSequenceInfo()) {
			if (prefixOnly) {
				if (_i.getHeader().startsWith(header)) info = _i;
			} else {
				if (_i.getHeader().equals(header)) info = _i;
			}
		}	
		StringBuffer sequence = new StringBuffer();
		if (info == null) return sequence;
        fastaFile.seek(info.getOffset());
        String line = fastaFile.readLine();
        if (prefixOnly) {
        	if (!line.startsWith(">"+header)) throw new InvalidIndexException("Sequence "+info.getSeqId()+" not at offset "+info.getOffset(),index);     
        } else {
        	if (!line.equals(">"+header)) throw new InvalidIndexException("Sequence "+info.getSeqId()+" not at offset "+ info.getOffset(),index);
        }
		while ((line = fastaFile.readLine()) != null) {
            if (line.startsWith(">"))
            	break;
            else
                sequence.append(line.trim());
        }
        return sequence;
    }
	
	public BaseStream getSequenceAsStream(String id) throws IOException, InvalidIndexException, ParseException {
		checkAndUpdateIndex();
		FastaDatabaseIndex.SequenceInfo info = null;
		for (SequenceInfo _i : index.getSequenceInfo()) {
			String prefix = _i.getHeader().split("\\s+")[0];
			if (prefix.equals(id)) {
				info = _i;
				break;
			}
		}
		if (info == null) return null;
		final BufferedRandomAccessFile file = new BufferedRandomAccessFile(fileFastaDatabase, "r");
		file.seek(info.getOffset());
		String line = file.readLine();
		if (!line.startsWith(">"+id)) throw new InvalidIndexException("Sequence "+info.getSeqId()+" not at offset "+info.getOffset(),index);     
		return new MyBaseStream(file); 
	}
	
	
	@Override
	public HashSet<String> writeSequencesByHeaders(File fileOutputFileFasta, HashSet<String> setOfHeadersWithoutGTSign, boolean keepNewlineFormatting, boolean truncateAfterFirstSpace) throws IOException {
		return super.writeSequencesByHeaders(fileOutputFileFasta, setOfHeadersWithoutGTSign, keepNewlineFormatting, truncateAfterFirstSpace);
    }
	
	@Override
	public void generatePaddedVersionToFile(int lengthOfPadOnFivePrimeEnd, int lengthOfPadOnThreePrimeEnd, char characterToPadWith, File fileOutputFastaPadded) throws IOException {
		super.generatePaddedVersionToFile(lengthOfPadOnFivePrimeEnd, lengthOfPadOnThreePrimeEnd, characterToPadWith, fileOutputFastaPadded);
    }
	
	@Override
	public void generateTruncatedVersionToFile(int lengthOfTruncatedReads,
            File fileOutputFasta) throws IOException {
		super.generateTruncatedVersionToFile(lengthOfTruncatedReads, fileOutputFasta);
	}

	@Override
    public void generateSplitVersionsToFiles(int lengthOfFirstPart, int lengthOfSecondPart,
            String stringToAppendToStartOfSecondPart,
            File fileOutputFastaFirstPart, File fileOutputFastaSecondPart) throws IOException {
    	super.generateSplitVersionsToFiles(lengthOfFirstPart, lengthOfSecondPart, stringToAppendToStartOfSecondPart, fileOutputFastaFirstPart, fileOutputFastaSecondPart);
	}
	
	@Override
    public ArrayList<HashSet<String>> getSetsOfSequenceHeadersEachTotalingLessThanNMegabases(long maxRefSequenceLength) throws Exception {
		return super.getSetsOfSequenceHeadersEachTotalingLessThanNMegabases(maxRefSequenceLength);
    }
	
	@Override
	public SortedMap<String, Long> getSequenceLengths(long max)
			throws Exception {
		if (max < 0) return getSequenceLengths();
		SortedMap<String, Long> map = getSequenceLengths();
		while (map.size() > max) {
			map.remove(map.firstKey());
		}
		return map;
	}
	
	@Override
    public SortedMap<String, Long> getSequenceLengths() throws IOException, ParseException, InvalidIndexException {
		checkAndUpdateIndex();
        TreeMap<String, Long> mapHeaderToSequenceLength = new TreeMap<String, Long>();
        for (SequenceInfo info : index.getSequenceInfo())
        	mapHeaderToSequenceLength.put(info.getHeader(), info.getLength());
        
        return mapHeaderToSequenceLength;
    }
	
	@Override
    public long getTotalLengthOfSequence() throws IOException, ParseException, InvalidIndexException {
		checkAndUpdateIndex();
        return index.getBaseCount();
    }
	
	public static void test(String[] args) throws Exception {
		File fastaFile = new File(args[0]);
		boolean doMemoryIntensiveTests = args.length > 1 && Boolean.parseBoolean(args[1]);
		File indexFile = new File("myindex");
		
		FastaDatabase database1 = new FastaDatabase(fastaFile);
		FastaDatabase database2 = new IndexedFastaDatabase(fastaFile, indexFile);
		System.out.println("getNumberOfSequencesInDatabase()");
		//if (database1.getNumberOfSequencesInDatabase() != database2.getNumberOfSequencesInDatabase() )
		//	System.out.println("getNumberOfSequencesInDatabase() test failed: "+database1.getNumberOfSequencesInDatabase() + " " + database2.getNumberOfSequencesInDatabase());

		if (doMemoryIntensiveTests) {
			System.out.print("nextSequence()");
			for (int i=0; i<database2.getNumberOfSequencesInDatabase(); i++) {
				System.out.print(".");
				SequenceHeaderPair pair1 = database1.nextSequence(false);	
				SequenceHeaderPair pair2 = database2.nextSequence(false);
				if (!pair1.getHeader().equals(pair2.getHeader())) 
					System.out.println("\nnextSequence() test failed:  '"+pair1.getHeader() +"' != '"+ pair2.getHeader() + "'");
				if (!pair1.getSequence().equals(pair2.getSequence()))
					System.out.println("\nnextSequence() test failed:  '"+pair1.getHeader() +"' != '"+ pair2.getHeader() + "'");
			}
			System.out.println();
		}
		
		if (doMemoryIntensiveTests) {
			System.out.println("nextSequence() with reset()");
			database1.reset();
			database2.reset();
			for (int i=0; i<2; i++) {
				if (i == 1) {
					database1.reset();
					database2.reset();
				}
				SequenceHeaderPair pair1 = database1.nextSequence(false);
				SequenceHeaderPair pair2 = database2.nextSequence(false);
				if (!pair1.getHeader().equals(pair2.getHeader())) 
					System.out.println("nextSequence()/reset() test failed:  '"+pair1.getHeader() +"' != '"+ pair2.getHeader() + "'");
				if (!pair1.getSequence().equals(pair2.getSequence()))
					System.out.println("nextSequence()/reset() test failed:  '"+pair1.getHeader() +"' != '"+ pair2.getHeader() + "'");
			}
		}
		
		System.out.println("getListOfHeaders()");
		List<String> headers1 = database1.getListOfHeaders(false, false);
		List<String> headers2 = database2.getListOfHeaders(false, false);
		if (!headers1.equals(headers2))
			System.out.println("getListOfHeaders() test failed");
	
		if (doMemoryIntensiveTests) {
			if (database2.getTotalLengthOfSequence() < 1000000000) {
				System.out.println("getHeaderToSequenceMap()");
				Map<String, String> map1 = database1.getHeaderToSequenceMap(false);
				Map<String, String> map2 = database2.getHeaderToSequenceMap(false);
				if (!map1.equals(map2))
					System.out.println("getHeaderToSequenceMap() test failed");
			} else {
				System.out.println("Skipping getHeaderToSequenceMap().  Database too big.");
			}
		}
		
		if (doMemoryIntensiveTests) {
			System.out.println("getSequenceByHeader()");
			int count = 0;
			List<String> headerCopy = headers1;
			if (headerCopy.size() > 100) {
				headerCopy = new ArrayList<String>(headers1);
				Collections.shuffle(headerCopy);
				headerCopy = headerCopy.subList(0, 80);
				headerCopy.addAll(headers1.subList(0,10));
				headerCopy.addAll(headers1.subList(headers1.size()-11, headers1.size() - 1));
			}
			for (String header : headerCopy) {
				String seq1 = database1.getSequenceByHeader(header).toString();
				String seq2 = database2.getSequenceByHeader(header).toString();
				if (count > 0) System.out.print("\b\b\b\b\b\b");
				System.out.printf("%6d", ++count);
				if (!seq1.equals(seq2)) System.out.println("getSequenceByHeader() test failed: "+header+":"+seq1.length()+" "+seq2.length());
			}
			System.out.println();
		}
		
		if (doMemoryIntensiveTests) {
			System.out.println("getSequenceByHeaderPrefix()");
			int count = 0;
			List<String> headerCopy = headers1;
			if (headerCopy.size() > 100) {
				headerCopy = new ArrayList<String>(headers1);
				Collections.shuffle(headerCopy);
				headerCopy = headerCopy.subList(0, 80);
				headerCopy.addAll(headers1.subList(0,10));
				headerCopy.addAll(headers1.subList(headers1.size()-11, headers1.size() - 1));
			}
			for (String header : headerCopy) {
				String prefix = header.split("\\s+")[0];
				String seq1 = database1.getSequenceByHeaderPrefix(prefix).toString();
				String seq2 = database2.getSequenceByHeaderPrefix(prefix).toString();
				if (count > 0) System.out.print("\b\b\b\b\b\b");
				System.out.printf("%6d", ++count);
				if (!seq1.equals(seq2)) System.out.println("getSequenceByHeaderPrefix() test failed: "+prefix+":"+seq1.length()+" "+seq2.length());
			}
			System.out.println();
		}
		
		System.out.println("getSequenceLengths()");
		SortedMap<String, Long> smap1 = database1.getSequenceLengths();
		SortedMap<String, Long> smap2 = database2.getSequenceLengths();
		
		if (!smap1.keySet().equals(smap2.keySet())) System.out.println("getSequenceLengths() keys test failed. "+smap1.size()+" "+smap2.size());
		if (!smap1.equals(smap2)) System.out.println("getSequenceLengths() test failed. "+smap1.size()+" "+smap2.size());
		
		System.out.println("getTotalLengthOfSequence()");
		if (database1.getTotalLengthOfSequence() != database2.getTotalLengthOfSequence())
			System.out.println("getTotalLengthOfSequence() test failed");
	}
	
	public static void main(String[] args) {
		try {
			test(args);
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
}		

class MyBaseStream implements BaseStream {

	private int baseCount = 0;
	private Character nextBase = null;
	private boolean hasMore = true;
	private BufferedRandomAccessFile file;

	public MyBaseStream(BufferedRandomAccessFile file) {
		this.file = file;
	}

	@Override
	public Character nextBase() throws IOException{
		if (!hasMore) return null;

		Character baseToReturn = null;

		while (hasMore && (nextBase == null || baseToReturn == null)) {
			if (baseToReturn == null) {
				baseToReturn = nextBase;
				nextBase = null;
			}
			int byt = file.read();
			if (byt < 0 || (char)byt == '>')
				hasMore = false;
			else if (!Character.isWhitespace(byt))
				nextBase = (char)byt;
		}
		if (baseToReturn != null) baseCount++;
		return baseToReturn;
	}

	@Override
	public int getBaseCount() {
		return baseCount;
	}

	@Override
	public boolean hasMoreBases() {
		return hasMore;
	}

	@Override
	public void close() throws IOException {
		file.close();
	}
}
