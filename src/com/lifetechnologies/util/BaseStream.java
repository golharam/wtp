package com.lifetechnologies.util;

import java.io.IOException;

/**
 * Represents a stream of base chars.  Keeps a tally of bases iterated.
 * @author mullermw
 *
 */
public interface BaseStream {
	
	/**
	 * 
	 * @return true if this stream has more bases.
	 */
	public boolean hasMoreBases();
	
	/**
	 * 
	 * @return the next base in the stream.
	 * @throws IOException
	 */
	public Character nextBase() throws IOException ;
	
	
	/**
	 * 
	 * @return count of bases read from this stream.
	 */
	public int getBaseCount();
	
	/**
	 * close the stream
	 * @throws IOException
	 */
	public void close() throws IOException;
	
	static class EmptyBaseStream implements BaseStream {
		@Override
		public boolean hasMoreBases() {
			return false;
		}
		
		@Override
		public void close() throws IOException {}
		
		@Override
		public int getBaseCount() {
			return 0;
		}
		
		@Override
		public Character nextBase() throws IOException {
			return null;
		}
	}
}
