package com.lifetechnologies.solid.wt.mapper;

/**
 * A Task in the WT mapping pipeline.
 * @author mullermw
 *
 */
public interface Task {

	/**
	 * Has this task been completed successfully?  Must not depend on
	 * prior execution of doTask()
	 * @return
	 * @throws TaskException
	 */
	public boolean isDone() throws TaskException;
	
	/**
	 * Perform this task.
	 * @throws TaskException
	 */
	public void doTask() throws TaskException;
	
	/**
	 * Whether this task should be executed remotely.
	 */
	public void setRemote(boolean remote);
	
	/**
	 * 
	 * @return true if this is a remote task.
	 */
	public boolean isRemote();
	
}
