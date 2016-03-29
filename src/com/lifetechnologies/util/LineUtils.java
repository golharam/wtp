package com.lifetechnologies.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.text.Format;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import java.text.ParseException;

public class LineUtils {

	/**
	 * Sort lines and write to a file.
	 * @param in
	 * @param out
	 * @param lineComparator
	 * @throws IOException
	 */
	public static void sortLines(Reader in, File out, Comparator<String> lineComparator) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		File tmp = null;
		try {
			if (in == null) return;
			if (out == null) return;
			List<String> headerLines = new ArrayList<String>();
			List<String> lines = new ArrayList<String>();
			
			reader = new BufferedReader(in);
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().startsWith("#")) {
					headerLines.add(line);
				} else {
					lines.add(line);
				}
			}
			
			if (lineComparator == null) {
				Collections.sort(lines);
			} else {
				Collections.sort(lines, lineComparator);
			}
			
			tmp = File.createTempFile(out.getName(), "tmp", out.getParentFile());
			tmp.deleteOnExit();
			
			writer = new PrintWriter(new BufferedWriter(new FileWriter(tmp)));
			for (String line : headerLines)
				writer.println(line);
			for (String line : lines) 
				writer.println(line);
			
			tmp.renameTo(out);
		} finally {
			if (writer != null) writer.close();
			if (tmp != null && tmp.exists()) tmp.delete();
		}
	}
	
	/**
	 * Sort lines by first parsing as objects, 
	 * using format to parse the line as an object 
	 * and format the object as a line.
	 * if lineComparator is null objects are sorted by their natual ordering.
	 * @param <T>
	 * @param in
	 * @param out
	 * @param lineComparator
	 * @param format
	 * @throws IOException
	 * @throws ParseException
	 */
	@SuppressWarnings("unchecked")
	public static <T extends Comparable> void sortLines(Reader in, File out, Comparator<T> lineComparator, Format format) throws IOException, ParseException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		File tmp = null;
		try {
			if (in == null) return;
			if (out == null) return;
			List<String> headerLines = new ArrayList<String>();
			List<T> objs = new ArrayList<T>();
			
			reader = new BufferedReader(in);
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().startsWith("#")) {
					headerLines.add(line);
				} else {
					objs.add((T)format.parseObject(line));
				}
			}
			
			if (lineComparator == null) {
				Collections.sort(objs);
			} else {
				Collections.sort(objs, lineComparator);
			}
			
			tmp = File.createTempFile(out.getName(), "tmp", out.getParentFile());
			tmp.deleteOnExit();
			
			writer = new PrintWriter(new BufferedWriter(new FileWriter(tmp)));
			for (String line : headerLines)
				writer.println(line);
			for (T obj : objs) 
				writer.println(format.format(obj));
			
			tmp.renameTo(out);
		} finally {
			if (writer != null) writer.close();
			if (tmp != null && tmp.exists()) tmp.delete();
		}
	}
	
	public static void sortLines(Reader in, File out, Format format) throws IOException, ParseException {
		sortLines(in, out, null, format);
	}
	
	/**
	 * Merge lines from several readers.  Lines in these readers should already be sorted.
	 * @param src
	 * @param dest
	 * @param comp
	 * @throws IOException
	 */
	public static void mergeLines(Collection<Reader> src, PrintWriter dest, Comparator<String> comp) throws IOException {
		
		if (dest == null) return;
		if (src == null) return;
		if (src.isEmpty()) return; 
		
		/*
		 * line2reader keeps the lines (keys) sorted.  The values are the readers from which the lines
		 * were read.
		 */
		SortedMap<String, BufferedReader> line2reader =  new TreeMap<String, BufferedReader>(comp);
		
		//Initialize line2reader
		for (Reader reader : src) {
			BufferedReader bufferedReader = reader instanceof BufferedReader ? (BufferedReader)reader : new BufferedReader(reader);
			String line = bufferedReader.readLine();
			if (line != null)
				line2reader.put(line, bufferedReader);
		}
		int lineCount = 0;
		long start = System.currentTimeMillis();
		while (line2reader.isEmpty() == false) {
			//Print the lowest line.
			String line = line2reader.firstKey();
			dest.println(line);
			//Get and remove the corresponding reader.
			BufferedReader reader = line2reader.get(line);
			line2reader.remove(line);
			//Read the next line from the reader.
			line = reader.readLine();
			//map the line to it's reader.
			if (line != null) line2reader.put(line, reader);
			if (++lineCount % 10000 == 0) {
				long now = System.currentTimeMillis();
				System.out.printf("Printing %d lines took %d millis.\n", lineCount, now-start);
			}
		}
	}
	
	/**
	 * Merge sorted lines from several readers by first sorting as objects 
	 * using format to parse line as an object and format object as a line.
	 * @param <T>
	 * @param src
	 * @param dest
	 * @param comp
	 * @param format
	 * @throws IOException
	 */
	@SuppressWarnings("unchecked")
	public static <T> void mergeLines(Collection<Reader> src, PrintWriter dest, Comparator<T> comp, Format format) throws IOException, ParseException {
		
		SortedMap <T, BufferedReader> obj2reader = comp == null ? new TreeMap<T, BufferedReader>() : new TreeMap<T, BufferedReader>(comp);
		//Initialize line2reader
		for (Reader reader : src) {
			BufferedReader bufferedReader = reader instanceof BufferedReader ? (BufferedReader)reader : new BufferedReader(reader);
			String line = bufferedReader.readLine();
			if (line != null)
				obj2reader.put((T)format.parseObject(line), bufferedReader);
		}
		
		while (obj2reader.isEmpty() == false) {
			//Print the lowest line.
			T obj = obj2reader.firstKey();
			dest.println(format.format(obj));
			//Get and remove the corresponding reader.
			BufferedReader reader = obj2reader.get(obj);
			obj2reader.remove(obj);
			//Read the next line from the reader.
			String line = reader.readLine();
			//map the line to it's reader.
			if (line != null)
				obj2reader.put((T)format.parseObject(line), reader);
		}
	}
	
	public static void mergeLines(Collection<Reader> src, PrintWriter dest, Format format) throws IOException, ParseException {
		mergeLines(src, dest, null, format);
	}
}

