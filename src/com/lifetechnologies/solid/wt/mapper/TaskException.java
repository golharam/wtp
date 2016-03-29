package com.lifetechnologies.solid.wt.mapper;

public class TaskException extends Exception {
	
	public static final long serialVersionUID = 1;
	
	public TaskException() {
		super();
	}
	
	public TaskException(String message) {
		super(message);
	}
	
	public TaskException(Throwable cause) {
		super(cause);
	}
	
	public TaskException(String message, Throwable cause) {
		super(message, cause);
	}
	

}
