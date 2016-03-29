package com.lifetechnologies.util;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

/**
 * Basic operations on Nucleotide sequence characters, including IUB codes.
 * @author mullermw
 *
 */
public class SeqUtils {

	private static String[] basesTrStrings = { "ACTGUMRWSYKVHDBXNactgumrwsykvhdbxn.", "TGACAKYWSRMBDHVXXtgacakywsrmbdhvxx."};
	public static final Map<Character, Character> baseTranslation;
	
	static {
		HashMap<Character, Character> tr = new HashMap<Character, Character>();
		char[] fromArr = basesTrStrings[0].toCharArray();
		char[] toArr   = basesTrStrings[1].toCharArray();
		for (int i=0; i<fromArr.length; i++) {
			tr.put(fromArr[i], toArr[i]);
		}
		baseTranslation = Collections.unmodifiableMap(tr);
	}
	
	private SeqUtils() {}
	
	/**
	 * 
	 * @param c
	 * @return return the base complement of c.
	 */
	public static final char complement(char c) {
		if (baseTranslation.containsKey(c)) return baseTranslation.get(c);
		return c;
	}
	
	/**
	 * 
	 * @param sequence 
	 * @return return the base complement sequence
	 */
	public static final String complement(String sequence) {
		if (sequence == null) return "";
		StringBuffer buffer = new StringBuffer();
		for (Character c : sequence.toCharArray())
			buffer.append(complement(c));
		return buffer.toString();
	}
	
	/**
	 * 
	 * @param sequence
	 * @return the reverse complement of sequence.
	 */
	public static final String reverseComplement(String sequence) {
		if (sequence == null) return "";
		return TextUtils.reverse(complement(sequence));

	}
	
	/**
	 * 
	 * @param buffer
	 * @return return the reverseComplement of buffer
	 */
	public static final StringBuffer reverseComplement(StringBuffer buffer) {
		if (buffer == null) return buffer;
		for (int i=0; i<buffer.length(); i++) {
			buffer.setCharAt(i, complement(buffer.charAt(i)));
		}
		buffer.reverse();
		return buffer;
	}
	
	public static final Pattern readPattern = Pattern.compile("^([A-Z])([0123\\.]{25,})$");
	
	/**
	 * Determine the read length in a mono-length csfasta file.
	 * @param csfasta
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 */
	public static int readLengthInCsfasta(InputStream csfasta) throws IOException, ParseException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(csfasta));
		int len = -1;
		int readCount = 0;
		int maxLines = 20;
		for (String line = reader.readLine(); line != null; ) {
			if (!line.startsWith(">")) continue;
			if (readCount++ > maxLines) break;
			line = reader.readLine();
			if (line == null) throw new ParseException("This stream is not csfasta.", readCount);
			Matcher matcher = readPattern.matcher(line);
			if (!matcher.matches()) throw new ParseException("This stream is not csfasta.", readCount);
			String colors = matcher.group(2);
			if (len < 0) len = colors.length();
			else if (len != colors.length()) throw new ParseException("This csfasta has reads of varying length.", readCount);
			line = reader.readLine();
			if (!line.startsWith(">")) throw new ParseException("This stream is not csfasta.", readCount);
		}
		return len;
	}
	
	@Test
	public void testReadLengthInCsfasta() throws IOException, ParseException {
		InputStream s = null;
		try {
			s = new FileInputStream("test/input/HBR.chr17_6.100k.mixed.csfasta");
			assertEquals(50, readLengthInCsfasta(s));
		} finally {
			s.close();
		}
	}
	
	public static void main(String[] args) {
		System.out.println(reverseComplement(basesTrStrings[0]));
		System.out.println(TextUtils.reverse(basesTrStrings[1]));
	}
}
