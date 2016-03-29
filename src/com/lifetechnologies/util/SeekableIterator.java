package com.lifetechnologies.util;

import java.util.Iterator;

/**
 * Extends java.util.Iterator with seekTo() and peek() functionality.
 * @author mullermw
 *
 * @param <T>
 */
public interface SeekableIterator<T> extends Iterator<T> {

	/**
	 * Jump to a value that matches prototype.
	 * The next call to next() will return the first value that is greater than
	 * or equal to prototype.
	 * @param prototype of the object we are seeking
	 */
	public void seekTo(T prototype);
	
	/**
	 * @return the next value but retain it for future calls to peek() or next()
	 */
	public T peek();
	
}
