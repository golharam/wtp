package com.lifetechnologies.solid.wt.cmd;

public class CmdRunner {

	private CmdRunner() {}
	
	/**
	 * This is the main method used by all the scripts in bin/
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		@SuppressWarnings("unchecked")
		Class<AbstractCmd> clazz = (Class<AbstractCmd>)Class.forName(args[0]);
		String[] newArgs = new String[args.length - 1];
		for (int i=1; i<args.length; i++) {
			newArgs[i-1] = args[i];
		}
		AbstractCmd cmd = clazz.newInstance();
		cmd.run(newArgs);
	}
	
}
