package com.lifetechnologies.solid.wt.mapper;

/**
 * Enumerates the ways to filter reads.
 * @author mullermw
 *
 */
public enum FilteringMode {
	/**
	 *  do not filter 
	 */ 
	OFF,
	
	/** 
	 * filter if one or more of the splits hits an entry in the filter reference 
	 */ 
	ONE_OR_MORE, 
	
	/**
	 * filter if both of the splits hit an entry in the filter reference.
	 */
	BOTH;
	
	public static FilteringMode parseFilteringMode(String s) {
		return FilteringMode.valueOf(s.toUpperCase());
	}
}
