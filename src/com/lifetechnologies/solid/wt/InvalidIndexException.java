package com.lifetechnologies.solid.wt;

/**
 * Thrown when and index is found to be inconsistent with it's subject.
 * @author mullermw
 *
 */
public class InvalidIndexException extends Exception {

	public static final long serialVersionUID = 1;
	private Object index;
	public String msg;
	
	public Object getIndex() {
		return index;
	}

	public void setIndex(Object index) {
		this.index = index;
	}

	/**
	 * 
	 * @param msg 
	 * @param index The index in question.
	 */
	public InvalidIndexException(String msg, Object index) {
		super(msg);
		this.index = index;
	}
}
