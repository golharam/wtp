package com.lifetechnologies.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.lifetechnologies.solid.wt.ExtendedReadMapping;
import com.lifetechnologies.solid.wt.ExtendedReadMappingByPositionComparator;
import com.lifetechnologies.solid.wt.SortedMaxFileIterator;

/**
 * An SeekableIterator that iterates over other (multiple) SeekableIterators.
 * @author mullermw
 *
 * @param <T>
 */
public class SeekablePolyIterator<T> implements SeekableIterator<T> {

	private Comparator<T> comparator;
	private List<SeekableIterator<T>> iterators;
	private SeekableIterator<T> nextIterator;
	
	/**
	 * @param comparator determines the ordering of values returned by next()
	 * @param iterators to be iterated over.
	 */
	public SeekablePolyIterator(Comparator<T> comparator, Set<SeekableIterator<T>> iterators) {
		this.comparator = comparator;
		this.iterators = new ArrayList<SeekableIterator<T>>(iterators);
		for (SeekableIterator<T> it : this.iterators)
			if (it.equals(this)) throw new IllegalArgumentException("Cannot iterate over self.");
	}
	
	@Override
	public boolean hasNext() {
		return readNextEntry(true) != null;
	}
	
	@Override
	public T next() {
		T value = readNextEntry(true);
		readNextEntry(false);
		return value;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("remove() not supported");
	}
	
	@Override
	public T peek() {
		return readNextEntry(true);
	}
	
	@Override
	public void seekTo(T prototype) {
		for (SeekableIterator<T> iterator : iterators)
			iterator.seekTo(prototype);
		return;
	}
	
	/*
	 * if reading from buffer, make sure the buffer is updated and return
	 * the value.
	 * if not reading from buffer, discard next value and buffer the next value
	 * return the next value.
	 */
	private T readNextEntry(boolean readFromBuffer) {
		
		//discard the buffer.
		if (readFromBuffer == false) {
			if (nextIterator != null) nextIterator.next();
			nextIterator = null;
		}
		
		if (nextIterator == null ) {	
			//Sort the iterators.
			Collections.sort(iterators, new Comparator<SeekableIterator<T>>() {
				@Override
				public int compare(SeekableIterator<T> o1, SeekableIterator<T> o2) {
					if (o1.equals(o2)) return 0;
					return comparator.compare(o1.peek(), o2.peek());
				}
			});
			//Buffer the iterator containing the first non-null value.
			for (SeekableIterator<T> iterator : iterators) {
				if (iterator.peek() != null) {
					nextIterator = iterator;
					break;
				}
			}
		}
		if (nextIterator == null) return null;
		return nextIterator.peek();
	}
	
	public static void main(String[] args) {
		try {
			File f1 = new File(args[0]);
			File f2 = new File(args[1]);
			ExtendedReadMapping prototype = SortedMaxFileIterator.createPrototype(Integer.parseInt(args[2]), Integer.parseInt(args[3]));
			Set<SeekableIterator<ExtendedReadMapping>> set = new HashSet<SeekableIterator<ExtendedReadMapping>>();
			set.add(new SortedMaxFileIterator(f1));
			set.add(new SortedMaxFileIterator(f2));
			SeekablePolyIterator<ExtendedReadMapping> it = new SeekablePolyIterator<ExtendedReadMapping>(new ExtendedReadMappingByPositionComparator(), set);
			it.seekTo(prototype);
			for (int i=0; i<20; i++) {
				System.out.println(it.next());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
