package com.lifetechnologies.util;

import java.text.ParseException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class MathUtils {
	
	private static Pattern byteStringPattern = Pattern.compile("^([\\d\\.]+)\\s*((?:k|m|g|t|p|e)?b?)$", Pattern.CASE_INSENSITIVE);
	public static final long KILOBYTE = 1024;
	public static final long MEGABYTE = 1048576;
	public static final long GIGABYTE = 1073741824;
	public static final long TERABYTE = 1099511627776L;
	public static final long PETABYTE = 1125899906842624L;
	public static final long EXABYTE =  1152921504606846976L;
	
	private MathUtils() {};
	
	/**
	 * Parse a String expressing a quantity of bytes (eg. 100mB, 5G, 100KB)
	 * to a number of bytes.
	 * @param s
	 * @return the number of bytes represented by the string.
	 * @throws ParseException
	 */
    public static long toBytes(String s) throws ParseException {
    	if (s == null) return 0;
    	s = s.trim();
    	if (s.matches("^\\d+$")) return Long.parseLong(s);
    	Matcher matcher = byteStringPattern.matcher(s);
    	double num;
    	if (matcher.matches()) {
    		num = Double.parseDouble(matcher.group(1));
    		String suffix = matcher.group(2);
    		
    		switch (suffix.toLowerCase().charAt(0)) {
    			
    			case 'b' :
    				//Leave it alone.
    				break;
    			case 'k' :
    				num *= KILOBYTE;
    				break;
    			case 'm' :
    				num *= MEGABYTE;
    				break;
    			case 'g' :
    				num *= GIGABYTE;
    				break;
    			case 't' :
    				num *= TERABYTE;
    				break;
    			case 'p' :
    				num *= PETABYTE;
    				break;
    			case 'e' :
    				num *= EXABYTE;
    			default :
    				throw new ParseException("Couldn't parse '"+s+"' as a byte string.", 0);
    		}
    	} else {
    		throw new ParseException("Couldn't parse '"+s+"' as a byte string.", 0);
    	}
    	return (long)num;
    }
    
}
