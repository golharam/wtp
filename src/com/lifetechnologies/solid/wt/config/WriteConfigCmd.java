package com.lifetechnologies.solid.wt.config;

import java.io.File;

public class WriteConfigCmd {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		ConfigWriter writer;
		if (args.length < 1)
			writer = new ConfigWriter(System.out);
		else 
			writer = new ConfigWriter(new File(args[0])); 
		try {
			for (ConfigKey key : ConfigKey.values())
				writer.writeDefaultConfig(key);
		} finally {
			writer.close();
		}
	}

}
