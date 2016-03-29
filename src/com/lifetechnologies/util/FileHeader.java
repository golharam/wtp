package com.lifetechnologies.util;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

/**
 * Composing, writing and reading File headers in WTP.
 * @author mullermw
 *
 */
public class FileHeader {

	public static final String HEADER_PREFIX = "#";
	public static final int MAX_LINES_WITHOUT_HEADER = 20;
	public static final String CREATED_BY_KEY = "created-by";
	public static final String SEQUENCE_SOURCE_FILE = "sequence-source-file";
	public static final String GENE_MODEL_SOURCE_FILE = "gene-model-source-file";
	public static final String SPLIT_KEY = "split";
	public static final String SPLIT_LENGTH_KEY = "split-length";
	public static final String MAX_MEMORY_KEY = "max-memory";
	public static final String JUNCTION_FLANK_SIZE = "junction-flank-size";
	public static final String JUNCTION_FASTA_FILE = "junction-fasta-file";
	public static final String FILE_SIZE = "file-size";
	public static final String JUNCTION_COUNT = "junction-count";
	
	/**
	 * Parse a header from an input stream.
	 * @param s
	 * @param headerPrefix The first character that appears in a header line (think '#')
	 * @return The header as a map.
	 * @throws IOException
	 */
	public static Map<String, Value> parseHeader(InputStream s, String headerPrefix) throws IOException {
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(s));
		StringBuffer header = new StringBuffer();
		int nonHeaderLineCount = 0;
		for (String line = reader.readLine(); line != null; line = reader.readLine()) {
			line = line.trim();
			
			if (!line.startsWith(headerPrefix)) {
				if (nonHeaderLineCount++ > MAX_LINES_WITHOUT_HEADER) break;
				else continue;
			} else { 
				nonHeaderLineCount = 0;
			}
			
			while (line.startsWith(headerPrefix)) {
				line = line.substring(1).trim();
			}
			header.append(line);
			header.append("\n");
		}
		Properties p = new Properties();
		p.load(new ByteArrayInputStream(header.toString().getBytes()));
		Map<String, Value> values = new HashMap<String, Value>();
		for (Map.Entry<Object, Object> entry: p.entrySet()) {
			String key = (String)entry.getKey();
			Value value = new FileHeader.Value(entry.getValue());
			values.put(key, value);
		}
		return values;
	}
	
	/**
	 * Convenience method.
	 * @param s
	 * @return
	 * @throws IOException
	 */
	public static Map<String, Value> parseHeader(InputStream s) throws IOException {
		return parseHeader(s, HEADER_PREFIX);
	}
	
	/**
	 * Parse a Java complieant properties stream.
	 */
	public static Map<String, Value> parseProperties(InputStream s) throws IOException {
		return parseHeader(s, "");
	}
	
	/**
	 * Convenience method.
	 * @param f
	 * @return
	 * @throws IOException
	 */
	public static Map<String, Value> parseHeader(File f) throws IOException {
		InputStream stream = null;
		try {
			stream = new FileInputStream(f);
			return parseHeader(stream);
		} finally {
			if (stream != null) stream.close();
		}
	}
	
	/**
	 * Convenience method.
	 * @param f
	 * @return
	 * @throws IOException
	 */
	public static Map<String, Value> parseProperties(File f) throws IOException {
		InputStream stream = null;
		try {
			stream = new FileInputStream(f);
			return parseProperties(stream);
		} finally {
			if (stream != null) stream.close();
		}
	}
	
	/**
	 * Write a header map to a stream.
	 * @param vals
	 * @param out
	 * @param comments
	 * @param commentPrefix
	 * @throws IOException
	 */
	public static void writeHeader(Map<?,?> vals, Writer out, String comments, String commentPrefix) throws IOException {
		if (vals == null) return;
		if (out == null) return;
		if (commentPrefix == null) commentPrefix = "#";
		Properties p = new Properties();
		if (vals instanceof Properties ) p = (Properties)vals;
		else
			//Copy values to p.
			for (Object obj : vals.entrySet()) {
				Map.Entry<?,?> entry = (Map.Entry<?,?>)obj;
				Object valueObj = entry.getValue();
				if (valueObj instanceof FileHeader.Value) {
					valueObj = ((FileHeader.Value)valueObj).toString();
				}
				if (valueObj == null) valueObj = "";
				p.put(entry.getKey(), valueObj);
			}
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		p.store(new PrintStream(baos), comments);
		BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(baos.toByteArray())));
		for (String line = reader.readLine(); line !=null; line = reader.readLine()) {
			out.write(commentPrefix+line+"\n");
		}
		out.flush();
	}
	
	/**
	 * convenience method
	 * @param vals
	 * @param out
	 * @param comments
	 * @param commentPrefix
	 * @throws IOException
	 */
	public static void writeHeader(Map<?,?> vals, OutputStream out, String comments, String commentPrefix) throws IOException {
		writeHeader(vals, new PrintWriter(out), comments);
	}
	
	/**
	 * convenience method
	 * @param vals
	 * @param out
	 * @param comments
	 * @throws IOException
	 */
	public static void writeHeader(Map<?,?> vals, OutputStream out, String comments) throws IOException {
		writeHeader(vals, out, comments, null);
	}
	
	/**
	 * convenience method
	 * @param vals
	 * @param out
	 * @param comments
	 * @throws IOException
	 */
	public static void writeHeader(Map<?,?> vals, Writer out, String comments) throws IOException {
		writeHeader(vals, out, comments, null);
	}

	/**
	 * A value in a header.  Convenient converstions to non-string datatypes.
	 * @author mullermw
	 *
	 */
	public static class Value {
		private String value;
		
		public Value(Object value) {
			this.value = value == null ? null : value.toString();
		}
				
		public boolean isNumeric() {
			if (value == null) return false;
			try {
				Double.parseDouble(value);
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		}
		public Integer intValue() {
			if (isNumeric()) return Integer.parseInt(value);
			return null;
		}
		
		public Long longValue() {
			if (isNumeric()) return Long.parseLong(value);
			return null;
		}
		
		public Double doubleValue() {
			if (isNumeric()) return Double.parseDouble(value); 
			return null;
		}
		
		@Override
		public String toString() {
			return value;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj == null) return false;
			if (this == obj) return true;
			if (!(obj instanceof FileHeader.Value)) return false;
			FileHeader.Value that = (FileHeader.Value)obj;
			if (this.toString() == null ) {
				if (that.toString() == null) return true;
				else return false;
			}
			if (that.toString() == null) return false;
			return this.toString().equals(that.toString());
		}
	}
}

