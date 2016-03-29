package com.lifetechnologies.util;

import java.io.ByteArrayOutputStream;
import java.lang.management.ManagementFactory;

import com.lifetechnologies.solid.wt.Utilities;

public class ProcessId {

	private static String processId = null;
	
	/**
	 * Get the process id if possible.
	 * @return the process id.
	 * @throws UnsupportedOperationException if process id cannot be retrieved on this platform.
	 */
    public static String getProcessId() throws UnsupportedOperationException {
    	
    	if (processId != null) return processId;
    	
    	/** This is a hack I found on the internet.  Currently works on Linux and Windows JVMs
    	 * Not guaranteed to work on all JVMs though so I double check with a pattern. **/
    	String runtimeMXBeanName = ManagementFactory.getRuntimeMXBean().getName();
    	if (!Utilities.isBlank(runtimeMXBeanName)) {
    		if (runtimeMXBeanName.matches("\\d+@\\S+")) {
    			processId = ManagementFactory.getRuntimeMXBean().getName().split("@")[0];
    		}
    	}
    	
    	if (processId != null) return processId;

    	String osName = System.getProperty("os.name");
    	/** Backup plan on non-windows:   use echo $PPID 
    	 * 
    	 * Even won't work with cygwin installed: $PPID is always 1.
    	 * 
    	 * */
    	if (!osName.toLowerCase().contains("windows")) {
	    	try {
		    	Process p = Runtime.getRuntime().exec(new String[] {"sh", "-c", "echo $PPID"});
		    	ByteArrayOutputStream baos = new ByteArrayOutputStream();
		    	Utilities.pump(p.getInputStream(), baos, true, true);
		    	p.waitFor();
		    	if (p.exitValue() == 0) processId = baos.toString().trim();
	    	} catch (Exception e) { }
    	}
    	if (processId == null) throw new UnsupportedOperationException("getProcessId() is not supported on this platform: " + osName);
    	return processId;
    }
    
	public static void main(String[] args) throws Exception {
		System.out.println(ProcessId.getProcessId());
		System.out.println(ProcessId.getProcessId());
	}
}
