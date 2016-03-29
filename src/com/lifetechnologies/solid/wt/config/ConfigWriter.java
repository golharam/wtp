package com.lifetechnologies.solid.wt.config;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

/**
 * Writes Configuration keys to an OutputStream as they should appear in 
 * a config file.  Each key is printed with a commented & wrapped description
 * Each key is printed with it's example value as a value.
 * 
 * Example:
 * 
 * ## Define the scheduling environment, should be 'pbs', 'lsf' or 'sge'
 * queue.sys=pbs
 * 
 * @author mullermw
 *
 */
public class ConfigWriter {

	private PrintWriter writer;
	
	public ConfigWriter(PrintWriter pw) {
		writer = pw;
	}
	
	public ConfigWriter(OutputStream os) {
		writer = new PrintWriter(os);
	}
	
	public ConfigWriter(File f) throws IOException {
		writer = new PrintWriter(new BufferedWriter(new FileWriter(f)));
	}
	
	public void close() {
		writer.close();
	}
	
	public void writeDefaultConfig(ConfigKey key) {
		String formattedDescription = key.getDescription().replaceAll("(.{77}\\s)", "$1\n");
		formattedDescription = formattedDescription.replaceAll("^", "## " );
		formattedDescription = formattedDescription.replaceAll("\n", "\n## ");
		writer.println(formattedDescription);
		writer.println(key.getQualifiedName()+"="+key.getExampleValue());
		writer.println();
	}
	
}
