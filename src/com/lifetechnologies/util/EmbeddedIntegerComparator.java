package com.lifetechnologies.util;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Compares strings based on the first parseable integer contained therein.
 * @author mullermw
 *
 */
public class EmbeddedIntegerComparator implements Comparator<String> {

	public static final EmbeddedIntegerComparator INSTANCE = new EmbeddedIntegerComparator();
	private static final Pattern digitPattern = java.util.regex.Pattern.compile("(\\d+)");
	@Override
	public int compare(String o1, String o2) {
		if (o1 == o2) return 0;
		if (o1.equals(o2)) return 0;
		Matcher matcher1 = digitPattern.matcher(o1);
		Matcher matcher2 = digitPattern.matcher(o2);
		int int1 = Integer.MAX_VALUE;
		int int2 = Integer.MAX_VALUE;
		if (matcher1.find())
			int1 = Integer.parseInt(matcher1.group());
		if (matcher2.find())
			int2 = Integer.parseInt(matcher2.group());
		if (int1 == int2) return o1.compareTo(o2);
		return int1 - int2;
	}
}
