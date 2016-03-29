package com.lifetechnologies.solid.wt;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReadMask {

	/**
	 * Represents a non-masked position in a ReadMask.
	 */
	public static final char INCLUDE_TOKEN = '1';
	
	/**
	 * Represents a masked position in a ReadMask.
	 */
	public static final char EXCLUDE_TOKEN = '0';
	
	/**
	 * Array corresponding to potential masked position.  false = masked.
	 */
	private Boolean[] model;
	
	/**
	 * Whether positions beyond those explicitly set in the model are masked.  false = masked.
	 */
	private boolean infinityMask = true;
	
	/**
	 * 
	 * A ReadMask is specifies whether a position is to be considered at
	 * each position of a read.  This constructor accepts String arguments
	 * of the following forms:
	 * 
	 * Bit Mask : A String of zeroes and ones to specify the masked and unmasked
	 *            positions respectively.
	 *            
	 *            Example:
	 *            
	 *            ReadMask m = new ReadMask("111111111111111100000");
	 *            
	 * Numeric Mask : A Comma separated list of numeric ranges identifying the
	 *                masked positions.  
	 *                
	 *            Example:
	 *            
	 *          		   "45" : position 45 is masked.
	 *					   "1..10, 15, 40..50" : mask positions 1 to 10, 15 and 40-50 
	 *            		   "45.." : position 45 and higher are masked.
	 *            
	 * 
	 * @param mask a String representation of the ReadMask. if mask is null, 
	 *        no positions will be masked.
	 */
	public ReadMask(String mask) {
		if (mask == null || mask.trim().equals("")) {
			model = new Boolean[0];
			return;
		}
		mask = mask.trim();
		if (mask.matches("^[01]+$")) {
			model = new Boolean[mask.length()];
			for (int i=0; i<mask.length(); i++)
				model[i] = mask.charAt(i) == INCLUDE_TOKEN;
		} else {
			List<MaskSegment> segments = new ArrayList<MaskSegment>();
			int length = 0;
			for (String segStr : mask.split("\\s*,\\s*")) {
				if (segStr.trim().isEmpty()) continue;
				MaskSegment segment = new MaskSegment(segStr);
				segments.add(segment);
				if (segment.getEnd() > length ) length = segment.getEnd();
				if (segment.getStart() > length ) length = segment.getStart();
			}
			model = new Boolean[length];
			for (int i = 0; i< model.length; i++) model[i] = true;
			for (MaskSegment segment : segments) {
				int end = segment.getEnd();
				if (end < 0) {
					end = model.length;
					infinityMask = false;
				}
				for (int i=segment.getStart()-1; i < end; i++)
					model[i] = false;
			}
		}
	}
	
	/**
	 * Is a position masked?
	 * @param i  position 
	 * @return true if it's masked.  false otherwise.
	 * @throws IllegalArgumentException if i is less than 1.
	 */
	public Boolean isMasked(int i) {
		if (i < 1) throw new IllegalArgumentException("Mask index is 1-based and cannot be less than 1.");
		if (i > model.length ) return !infinityMask;
		return !model[i-1];
	}
	
	/**
	 * Is a position not masked?
	 * @param i position
	 * @return true if it's not masked.  false otherwise.
	 */
	public Boolean isNotMasked(int i) {
		return !isMasked(i);
	}
	
	/**
	 * 
	 * @param length
	 * @return a string of length specified by the argument.  The string will
	 * consist of 0's and 1's corresponding to masked and unmasked positions
	 * respectively.  Ex: "1111111111111111100000"
	 */
	public String asBitString(int length) {
		StringBuffer buf = new StringBuffer();
		if (model != null)
			for (int i=1; i<=length; i++) {
				if (this.isMasked(i)) buf.append(EXCLUDE_TOKEN);
				else buf.append(INCLUDE_TOKEN);
			}
		return buf.toString();
	}
	
	/**
	 * 
	 * @return a bit string corresponding to the explicitly set masked positions.
	 * @see asBitString(length)
	 */
	public String asBitString() {
		return asBitString(model.length);
	}
	
	public ReadMask subMask(int from, int to) {
		return new ReadMask(this.asBitString(to).substring(from, to));
	}
	
	public int numberOfUnmaskedPositions(int length) {
		String bitString = this.asBitString(length);
		int count = 0;
		for (char c : bitString.toCharArray())
			if (c == '1') count++;
		return count;
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName().toString()+":["+this.asBitString()+","+(infinityMask ? 1 : 0)+"]";
	}

}

/**
 * 
 * Represents a segment of masked positions.  Not to be used outside ReadMask.  Positions are 1-based.
 * @author mullermw
 *
 */
class MaskSegment {
	
	/**
	 * the lowest (1-based) position of the segment.
	 */
	private int start;

	/**
	 * the length of the segment.  negative value indicates undefined length.
	 */
	private int length;

	/**
	 * Regular expression matching a segment.
	 */
	private static Pattern segmentPattern = Pattern.compile("^(\\d+)\\s*(?:\\.\\.|-)\\s*(\\d+)?$");
	
	/**
	 * Construct a segment from a compliant string.
	 * @param s
	 * @throws IllegalArgumentException if an invalid string is passed.
	 */
	public MaskSegment(String s) {
		if (s == null) throw new IllegalArgumentException("mask cannot be null");
		s = s.trim();
		if (s.matches("^\\d+$")) {
			start = Integer.parseInt(s);
			length = 1;
			return;
		}
		Matcher matcher = segmentPattern.matcher(s);
		if (matcher.matches()) {
			start = Integer.parseInt(matcher.group(1));
			if (start < 1) throw new IllegalArgumentException("Mask start positions are 1-based and cannot be less than one: "+s);
			if (matcher.group(2) == null) {
				length = -1;
			} else {
				int end = Integer.parseInt(matcher.group(2));
				if (end < start) {
					int temp = start;
					start = end;
					end = temp;
				}
				length = end - start + 1;
			}
			return;
		}
		throw new IllegalArgumentException("Unrecognized Mask Segment pattern: " + s);
	}
	
	/**
	 * 
	 * @return The starting position (lowest 1-based position) of the segment.
	 */
	public int getStart() {
		return start;
	}

	/**
	 * 
	 * @return Length of the segment. A negative value indicates an undefined end.
	 */
	public int getLength() {
		return length;
	}
	
	/**
	 * 
	 * @return end of the segment. A negative value indicates an undefined end.
	 */
	public int getEnd() {
		if (length < 0) return -1;
		return start + length - 1;
	}
	
	@Override 
	public String toString() {
		return getStart() + "-" + getEnd() + ":" + getLength();
	}
}
