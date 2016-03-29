package com.lifetechnologies.util;

import java.util.ArrayList;
import java.util.List;

public class TextUtils {

	private TextUtils() {}
	
	/**
	 * Split a line of comma separated values.  Quoted text is not split.  
	 * Backslash escaped commas and quotes are interpreted literally.
	 * @param value comma separated line of text
	 * @return list of substrings.
	 */
	public static List<String> splitCSV(String value) {
		List<String> values = new ArrayList<String>();
		if (value == null) return values;
		StringBuffer currValue = new StringBuffer();
		
		//Marks the position of an opening quote in currValue;
		//Non null value means we are inside quotes.
		Integer openQuoteMark = null; 
		
		for (int i=0; i<value.length(); i++) {
			char c = value.charAt(i);
			
			// Backslash escape, take next char literally.
			if (openQuoteMark == null && c == '\\') {
				if (++i < value.length()) currValue.append(value.charAt(i));
				continue;
			}
			
			if (isQuote(c)) {
				if (openQuoteMark == null) {
					// mark as a quoted region
					currValue.append(c);
					openQuoteMark = currValue.length() - 1;
					continue;
				} else if (c == currValue.charAt(openQuoteMark)) {
					// unmark quoted region.
					currValue.deleteCharAt(openQuoteMark);
					openQuoteMark = null;
					continue;
				} 
			}
			
			//Unquoted comma, end of value.
			if (c == ',' && openQuoteMark == null ) {
				values.add(currValue.toString());
				currValue.setLength(0);
				continue;
			}
			currValue.append(c);
		}
		values.add(currValue.toString());
		return values;
	}
	
	/**
	 * 
	 * @param c
	 * @return true if c is a ' or a "
	 */
	public static boolean isQuote(char c) {
		return c == '\'' || c == '"';
	}

	/**
	 * @param length
	 * @return return a length long string of space chars.
	 */
	public static String getSpaces(int length) {
		StringBuffer buf = new StringBuffer();
		for (int i=0; i<Math.abs(length); i++) {
			buf.append(" ");
		}
		return buf.toString();
	}
	
	/**
	 * 
	 * @param s
	 * @return the reversal of s
	 */
	public static String reverse(String s) {
		return new StringBuffer(s).reverse().toString();
	}
	
	public static void main(String[] args) {
		String[] strings = {"string one \\' \\\\i',s, her'e;" };
		for (String string : strings ) {
			System.out.println(string+"->"+splitCSV(string)+" ("+splitCSV(string).size()+")");
		}
	}
}
