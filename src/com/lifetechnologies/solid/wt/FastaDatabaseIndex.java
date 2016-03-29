package com.lifetechnologies.solid.wt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DateFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Represents a persistant index on a Fasta File. 
 * 
 * Example of an index file.
 * 
 * Table columns are description line, byte offset of '>' and seq length.
 * 
 * #target_path=/share/reference/genomes/chromFa/human.fa
 *	#file_size=3142044949
 * #date_modified=Nov 1, 2007 4:58:44 PM PDT
 * #line_length=50
 * #base_count=3080436051
 * chr1    0       247249719
 * chr2    252194720       242951149
 * chr3    500004898       199501827
 * chr4    703496768       191273063
 * chr5    898595299       180857866
 * chr6    1083070329      170899992
 * chr7    1257388327      158821424
 * chr8    1419386186      146274826
 * chr9    1568586515      140273252
 * chr10   1711665239      135374737
 * chr11   1849747478      134452384
 * chr12   1986888917      132349534
 * chr13   2121885449      114142980
 * chr14   2238311296      106368585
 * chr15   2346807260      100338915
 * chr16   2449152961      88827254
 * chr17   2539756768      78774742
 * chr18   2620107012      76117153
 * chr19   2697746516      63811651
 * chr20   2762834408      62435964
 * chr21   2826519099      46944323
 * chr22   2874402316      49691432
 * chrX    2925087584      154913754
 * chrY    3083099620      57772954
 * chrM    3142028040      16571
 * 
 * @author mullermw
 *
 */
public class FastaDatabaseIndex {
		
	private File target;
	private Long filesize;
	private Date dateModified;
	private Integer lineLength;
	private Long baseCount;
	private final ArrayList<SequenceInfo> sequenceInfos = new ArrayList<SequenceInfo>();
	private final Map<String, SequenceInfo> header2info = new HashMap<String, SequenceInfo>();
	private static final Pattern keyValuePattern = Pattern.compile("^#(\\S+)=(.*)");
	private static final String keyValuePrintf = "#%s=%s\n";
	private static final Pattern seqOffsetPattern = Pattern.compile("^(.+)\\t(\\d+)\\t(\\d+)");
	private static final String seqOffsetPrintf = "%s\t%s\t%s\n";
	private static final String KEY_TARGET_PATH = "target_path";
	private static final String KEY_FILE_SIZE = "file_size";
	private static final String KEY_DATE_MODIFIED = "date_modified";
	private static final String KEY_LINE_LENGTH = "line_length";
	private static final String KEY_BASE_COUNT = "base_count";
	private static final DateFormat dateFormat = DateFormat.getDateTimeInstance(DateFormat.MEDIUM, DateFormat.FULL, Locale.US);
	private static final Logger logger = Logger.getLogger(FastaDatabaseIndex.class.toString());
	
	protected FastaDatabaseIndex() {}
	
	protected FastaDatabaseIndex(InputStream s) throws IOException, ParseException {
		this();
		readIndex(s);
	}
	
	/**
	 *  New index from a fasta file.
	 * @param f the file to index.
	 * @return the index
	 * @throws IOException if an IO problem occurs in reading the fasta.
	 * @throws ParseException  if a problem occurs in parsing the fasta.
	 */
	public static FastaDatabaseIndex newIndexFromFasta(File f) throws IOException, ParseException {
		FastaDatabaseIndex index = new FastaDatabaseIndex();
		index.buildIndex(f);
		return index;
	}
	
	/**
	 *  New index from an input stream matching the proper format.
	 * @param s the stream to parse
	 * @return the index
	 * @throws IOException if an IO problem occurs in reading the stream.
	 * @throws ParseException if a problem occurs in parsing the stream.
	 */
	public static FastaDatabaseIndex newIndexFromIndexStream(InputStream s) throws IOException, ParseException {
		return new FastaDatabaseIndex(s);
	}
	
	/**
	 * Clear all fields of the index.
	 */
	public void clear() {
		target = null;
		filesize = null;
		dateModified = null;
		lineLength = null;
		baseCount = null;
		sequenceInfos.clear();
		header2info.clear();
	}
	
	/**
	 * Parse an input stream to update this index.
	 * @param instream the stream to parse
	 * @throws IOException if an IO error occurs in reading the stream.
	 * @throws ParseException if a problem occurs in parsing the index.
	 */
	public void readIndex(InputStream instream) throws IOException, ParseException {
		clear();
		BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
		for (String line = reader.readLine(); line != null; line = reader.readLine()) {
			if (Utilities.isBlank(line)) continue;
			Matcher matcher = keyValuePattern.matcher(line);
			if (matcher.matches()) {
				String key = matcher.group(1);
				String value = matcher.group(2);
				if (key.equals(KEY_TARGET_PATH))
					target = new File(value);
				else if (key.equals(KEY_FILE_SIZE))
					filesize = Long.valueOf(value);
				else if (key.equals(KEY_DATE_MODIFIED))
					dateModified = dateFormat.parse(value);
				else if (key.equals(KEY_LINE_LENGTH))
					lineLength = Integer.valueOf(value);
				else if (key.equals(KEY_BASE_COUNT))
					baseCount = Long.valueOf(value);
				continue;
			}
			matcher = seqOffsetPattern.matcher(line);
			if (matcher.matches()) {
				SequenceInfo info = new SequenceInfo(matcher.group(1), Long.parseLong(matcher.group(2)), Long.parseLong(matcher.group(3)));
				this.addSequenceInfo(info);
				continue;
			}
			throw new ParseException("Invalid line: "+line, 0);
		}
	}
	
	/**
	 * Read and parse a fasta file to update this index.
	 * @param file the fasta file to read.
	 * @throws IOException if a problem occurs in reading the file.
	 * @throws ParseException if a problem occurs in parsing the file.
	 */
	public void buildIndex(File file) throws IOException, ParseException {
		clear();
		if (file == null) return;
		this.target = file.getAbsoluteFile();
		this.filesize = file.length();
		this.dateModified = new Date(file.lastModified());
		BufferedReader reader = new BufferedReader(new FileReader(file));
		try {
			long offset = 0;
			long linecount = 0;
			Integer lastLineLength = null;
			this.baseCount = 0L;
			SequenceInfo curr = null;
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				linecount++;
				if (line.startsWith(">")) {
					if (curr != null)
						this.addSequenceInfo(curr);
					String header = line.substring(1);
					curr = new SequenceInfo(header, offset);
					logger.info("Processing: "+curr.getSeqId());
					lastLineLength = null;
				} else {
					if (this.lineLength == null)
						this.lineLength = line.length();
					else if (lastLineLength != null && lastLineLength != this.lineLength )
						throw new ParseException("Different line length in "+file+" line: "+ linecount, (int)offset );
					lastLineLength = line.length();
					if (curr != null) curr.length += line.length();
					this.baseCount += line.length();
				}
				offset += line.length() + 1;
			}
			if (curr != null)
				this.addSequenceInfo(curr);
		} finally {
			reader.close();
		}
	}
	
	/**
	 * Rebuild this index using the same target fasta file.
	 * @throws IOException if a problem occurs in IO.
 	 * @throws ParseException if a problem occurs in parsing.
	 */
	public void rebuildIndex() throws IOException, ParseException {
		this.buildIndex(this.target);
	}
	
	/**
	 * Write this index to a file in the proper format.
	 * @param writer to write to.
	 * @throws IOException if a IO problem occurs.
	 */
	public void writeIndex(PrintWriter writer) throws IOException {
		writer.printf(keyValuePrintf, KEY_TARGET_PATH, target.getAbsolutePath());
		writer.printf(keyValuePrintf, KEY_FILE_SIZE, filesize);
		writer.printf(keyValuePrintf, KEY_DATE_MODIFIED, dateFormat.format(dateModified));
		writer.printf(keyValuePrintf, KEY_LINE_LENGTH, lineLength);
		writer.printf(keyValuePrintf, KEY_BASE_COUNT, baseCount);
		for (SequenceInfo info : sequenceInfos )
			writer.printf(seqOffsetPrintf, info.getHeader(), info.getOffset(), info.getLength());
		writer.flush();
	}
	
	/**
	 * Write this index to a file.
	 * @param f
	 * @throws IOException
	 */
	public void writeIndex(File f) throws IOException {
		PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(f)));
		try {
			this.writeIndex(writer);
		} finally {
			writer.close();
		}
	}
	
	/**
	 * @return true if the index fields match the target fasta file.  Sequence
	 *         entries are not checked.
	 * @throws IOException
	 */
	public boolean isAccurate() throws IOException {
		return isAccurate(0);
	}
	
	/**
	 * 
	 * @param numSeqsToCheck The number of randomly selected sequences to check for proper indexing.
	 * @return true if the index fields match the target fasta file.
	 * @throws IOException
	 */
	public boolean isAccurate(int numSeqsToCheck) throws IOException {
		if (!target.exists()) return false;
		if (!target.isFile()) return false;
		if (target.length() != filesize) return false;
		Date currDateModified = new Date(target.lastModified());
		if (!currDateModified.equals(dateModified)) return false;
		if (numSeqsToCheck > this.getNumberOfSequences()) numSeqsToCheck = this.getNumberOfSequences();
		if (numSeqsToCheck > 0) {
			int numChecked = 0;
			BufferedRandomAccessFile file = new BufferedRandomAccessFile(target, "r");
			try {
				//Check 25 randomly selected sequences sequences.
				for (Iterator<SequenceInfo> it=header2info.values().iterator(); it.hasNext() && numChecked < numSeqsToCheck; ) {
					SequenceInfo info = it.next();
					file.seek(info.getOffset());
					String line = file.readLine();
					if (!line.startsWith(">"+info.getHeader())) return false;
					numChecked++;
				}
			} finally {
				file.close();
			}
		} 
		return true;
	}
	
	/**
	 * true if all the fields of the index match.
	 */
	public boolean equals(Object thatObj) {
		if (thatObj == null) return false;
		if (this == thatObj) return true;
		if (!(thatObj instanceof FastaDatabaseIndex)) return false;
		FastaDatabaseIndex that = (FastaDatabaseIndex)thatObj;
		if (this.dateModified == null)
			if (that.dateModified != null) return false;
		else
			if (!this.dateModified.equals(that.dateModified)) return false;

		if (this.filesize == null)
			if (that.filesize != null) return false;
		else
			if (!this.filesize.equals(that.filesize)) return false;

		if (this.target == null)
			if (that.target != null) return false;
		else
			if (!this.target.equals(that.target)) return false;
		
		if (this.lineLength == null)
			if (that.lineLength != null) return false;
		else
			if (!this.lineLength.equals(that.lineLength)) return false;
		
		if (this.baseCount == null)
			if (that.baseCount != null) return false;
		else
			if (!this.baseCount.equals(that.baseCount)) return false;
		
		if (!this.header2info.equals(that.header2info)) return false;
		
		return true;
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		FastaDatabaseIndex index = new FastaDatabaseIndex();
		index.dateModified = this.dateModified;
		index.filesize = this.filesize;
		index.target = this.target;
		index.baseCount = this.baseCount;
		index.lineLength = this.lineLength;
		for (SequenceInfo info : this.sequenceInfos)
			index.addSequenceInfo(info);
		return index;
	}
	
	public String toString() {
		StringWriter stringWriter = new StringWriter();
		try {
			writeIndex(new PrintWriter(stringWriter));
		} catch (IOException e) {
			stringWriter.write("ERROR");
		}
		return stringWriter.toString();
	}
	
	protected void addSequenceInfo(SequenceInfo info) {
		sequenceInfos.add(info);
		header2info.put(info.getHeader(), info);
	}
	
	public File getTarget() {
		return target;
	}

	public Long getFilesize() {
		return filesize;
	}

	public Date getDateModified() {
		return dateModified;
	}

	public Integer getLineLength() {
		return lineLength;
	}

	public Long getBaseCount() {
		return baseCount;
	}
	
	public int getNumberOfSequences() {
		return sequenceInfos.size();
	}
	
	public SequenceInfo getSequenceInfo(int i) {
		return sequenceInfos.get(i);
	}

	public SequenceInfo getSequenceInfo(String header) {
		return header2info.get(header);
	}
	
	@SuppressWarnings("unchecked")
	public List<SequenceInfo> getSequenceInfo() {
		return (List<SequenceInfo>)sequenceInfos.clone();
	}

	public static void main(String[] args) {
		try {
			logger.setLevel(Level.INFO);
			FastaDatabaseIndex index = FastaDatabaseIndex.newIndexFromFasta(new File(args[0]));
			File indexFile = new File("myindex");
			PrintWriter writer = new PrintWriter(new FileOutputStream(indexFile));
			index.writeIndex(writer);
			System.out.println("index.isAccurate()="+index.isAccurate(25));
			
			FastaDatabaseIndex anotherIndex = FastaDatabaseIndex.newIndexFromIndexStream((new FileInputStream(indexFile)));
			System.out.println("index.equals(anotherIndex)="+index.equals(anotherIndex));
			System.out.println("anotherIndex.isAccurate(true)="+anotherIndex.isAccurate(25));
			
			File anotherIndexFile = new File("myindex.copy");
			writer = new PrintWriter(new FileOutputStream(anotherIndexFile));
			anotherIndex.writeIndex(writer);
			anotherIndex = FastaDatabaseIndex.newIndexFromIndexStream((new FileInputStream(anotherIndexFile)));
			System.out.println("anotherIndex.equals(index)="+anotherIndex.equals(index));
			System.out.println("anotherIndex.isAccurate(true)="+anotherIndex.isAccurate(25));
			
			FastaDatabaseIndex badIndex = (FastaDatabaseIndex)index.clone();
			badIndex.dateModified = new Date();
			System.out.println("!badIndex.isAccurate(true)="+!badIndex.isAccurate(25));
			badIndex.writeIndex(new PrintWriter(new FileOutputStream(new File("myindex.bad"))));

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Basic information about a Fasta Sequence.
	 * header is the description line of the record.
	 * offset is the byte offset of the ">" character in the record.
	 * length is the base length of the record.
	 * @author mullermw
	 *
	 */
	class SequenceInfo implements Comparable<SequenceInfo>{
		
		private String header;
		private Long offset = 0L;
		private Long length = 0L;
		
		
		public SequenceInfo(String header, Long offset) {
			super();
			this.header = header;
			this.offset = offset;
		}
		public SequenceInfo(String header, Long offset, Long length) {
			super();
			this.header = header;
			this.offset = offset;
			this.length = length;
		}
		public String getHeader() {
			return header;
		}
		public void setHeader(String header) {
			this.header = header;
		}
		public String getSeqId() {
			if (header == null) return null;
			return header.split("\\s+")[0];
		}
		public Long getOffset() {
			return offset;
		}
		public void setOffset(Long offset) {
			this.offset = offset;
		}
		public Long getLength() {
			return length;
		}
		public void setLength(Long length) {
			this.length = length;
		}
		
		public boolean equals(Object thatObj) {
			if (thatObj == null) return false;
			if (this == thatObj) return true;
			SequenceInfo that = (SequenceInfo)thatObj;
			if (this.header == null) {
				if (that.header != null) return false;
			} else {
				if (!this.header.equals(that.header)) return false;
			}
			if (this.offset == null) {
				if (that.offset != null) return false;
			} else {
				if (!this.offset.equals(that.offset)) return false;
			}
			if (this.length == null) {
				if (that.length != null) return false;
			} else {
				if (!this.length.equals(that.length)) return false;
			}
			return true;
		}
		
		@Override
		protected Object clone() throws CloneNotSupportedException {
			return new SequenceInfo(this.header, this.offset, this.length);
		}
		
		@Override
		public int compareTo(SequenceInfo that) {
			if (that == null) return 1;
			if (this.getOffset() > that.getOffset()) return 1;
			if (this.getOffset() < that.getOffset()) return -1;
			return 0;
		}
		
		@Override
		public String toString() {
			return Utilities.join(":", this.getClass(), this.header, this.offset, this.length );
		}
	}
}

