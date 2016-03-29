package com.lifetechnologies.solid.wt.cmd;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

import com.lifetechnologies.solid.wt.config.ComplexOption;
import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.solid.wt.config.ValidatingOption;
import com.lifetechnologies.solid.wt.config.WTConfigUtils;
import com.lifetechnologies.solid.wt.config.validator.MaxValuesValidator;
import com.lifetechnologies.solid.wt.config.validator.ReadableFileValidator;
import com.lifetechnologies.solid.wt.config.validator.ValidatorResult;
import com.lifetechnologies.util.TextUtils;

/**
 * Superclass for WT command line programs.
 * @author mullermw
 *
 */
public abstract class AbstractCmd {
	
	private String cmdName = "";
	private Options options;
	Map<ConfigKey, Set<String>> config = new HashMap<ConfigKey, Set<String>>();
	static final Logger logger = Logger.getLogger(AbstractCmd.class.toString());
	private static HelpFormatter helpFormatter = new MyHelpFormatter();
	private final String command;
	
	static {
		helpFormatter.setOptionComparator(new Comparator<Option>() {
			@Override
			public int compare(Option o1, Option o2) {
				if (o1 instanceof ComplexOption && o2 instanceof ComplexOption) {
					ComplexOption co1 = (ComplexOption)o1;
					ComplexOption co2 = (ComplexOption)o2;
					return co1.getKey().compareTo(co2.getKey());
				}
				return o1.getLongOpt().compareTo(o2.getLongOpt());
			}
		});
	}
	
	public AbstractCmd(String name) {
		if (name == null) name = "";
		this.cmdName = name;
		command = System.getProperty("command", this.getClass().getSimpleName());
	}
	
	public void preRun(String[] args) throws ParseException, IOException{
		System.setProperty("java.awt.headless", "true"); 
		LogManager logManager = LogManager.getLogManager();
		logManager.reset();
		Logger rootLogger = logManager.getLogger("");
		ConsoleHandler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.SEVERE);
		rootLogger.setLevel(Level.INFO);
		rootLogger.addHandler(consoleHandler);

		logger.info("Starting "+getCmdName());
		//Contains command line options
		
		logger.info("Processing command line options");
		Options options = this.getOptions();

		//Parse the command line.
		Parser parser = new PosixParser();
		CommandLine cmdLine = null;
		try {
			cmdLine = parser.parse(options, args);
		} catch(org.apache.commons.cli.ParseException e) {
			System.out.println(e.getMessage());
			helpFormatter.printUsage(new PrintWriter(System.out, true) , MyHelpFormatter.TERMINAL_WIDTH, command, options);
			System.out.printf("type '%s -h' for help\n", command);
			logger.severe("Error processing command line options");
			System.exit(1);
		}
		
		if (cmdLine.hasOption("help")) {
			helpFormatter.printHelp(command, options);
			System.exit(0);
		}
		
		//Initialize the config with defaults.
		Map<ConfigKey, Set<String>> config = WTConfigUtils.newConfig();
		for (Object obj : options.getOptions()) {
			if (obj instanceof ComplexOption == false) continue;
			ComplexOption option = (ComplexOption)obj;
			ConfigKey key = option.getKey();
			if (key == null) {
				System.out.println(option);
				System.exit(1);
			}
			if (key.hasDefault()) config.put(key, WTConfigUtils.asSingleton(key.getDefaultValue()));
		}

		//Parse config files.
		logger.info("Processing config files");
		if (cmdLine.hasOption("config-file")) {
			ValidatingOption option = (ValidatingOption)options.getOption("config-file");
			ValidatorResult result = option.validateAll(cmdLine.getOptionValues("config-file"), config);
			if (result.isValid()) {
				File file = new File(cmdLine.getOptionValue("config-file"));
				WTConfigUtils.mapConfig(config, file);
			} else {
				StringBuffer bigMessage = new StringBuffer("config-file: \n");
				for (String message : result.getMessages()) {
					bigMessage.append(message.concat("\n"));
				}
				logger.severe(bigMessage.toString());
				System.exit(1);
			}
		} else {
			WTConfigUtils.mapConfig(config);
		}
		logger.fine(config.toString());
		
		//Remove non-relevant values from config.
		Set<ConfigKey> keysInUse = new HashSet<ConfigKey>();
		for (Object obj : options.getOptions()) {
			if (obj == null || !(obj instanceof ComplexOption)) continue;
			ComplexOption option = (ComplexOption)obj;
			keysInUse.add(option.getKey());
		}
		for (Iterator<ConfigKey> it = config.keySet().iterator(); it.hasNext();)
			if (!keysInUse.contains(it.next()))	it.remove();
		
		
		//Merge command-line options into config.
		for (Option option : cmdLine.getOptions()) {
			if (option instanceof ComplexOption) {
				ComplexOption defOpt = (ComplexOption)option;
				ConfigKey key = defOpt.getKey();
				if (key == null) continue;
				if (defOpt.hasArg()) {
					String[] valueArr = cmdLine.getOptionValues(option.getLongOpt());
					if (valueArr == null) valueArr = new String[0];
					config.put(key, new HashSet<String>(Arrays.asList(valueArr)));
				} else {
					config.put(key, new HashSet<String>(Arrays.asList(new String[] {"true"})));
				}
			}
		}
		
		//Dump config and exit.
		if (cmdLine.hasOption("print-parameters")) {
			List<ConfigKey> configKeys = new ArrayList<ConfigKey>(config.keySet());
			Collections.sort(configKeys);
			for (ConfigKey key : configKeys) {
				StringBuffer valueStr = new StringBuffer();
				for (String value : config.get(key))
					valueStr.append(",").append(value);
				valueStr.deleteCharAt(0);
				System.out.printf("%s='%s'\n", key.getQualifiedName(), valueStr);
			}
			System.exit(0);
		}
		
		//Validate config.
		logger.info("Validating config");
		boolean valid = true;
		for (Object optionObj : options.getOptions()) {
			Option option = (Option)optionObj;
			if (option instanceof ComplexOption) {
				ComplexOption complexOption = (ComplexOption)option;
				ConfigKey key = complexOption.getKey();
				if (key == null) continue;
				Set<String> values = config.get(key);
				ValidatorResult result = complexOption.validateAll(values, config);
				if (!result.isValid()) {
					valid = false;
					System.out.printf("  -%s  -%s  %s\n", key.asCommandLineShortOpt(), key.asCommandLineLongOpt(), key.getQualifiedName());
					for (String message : result.getMessages()) 
						System.err.printf("    %s\n", message);
					System.out.println();
				}
			}
		}

		if (!valid) {
			System.err.println("Invalid command line options");
			helpFormatter.printUsage(new PrintWriter(System.out, true) , MyHelpFormatter.TERMINAL_WIDTH, command, options);
			System.out.printf("type '%s -h' for help\n", command);
			System.exit(1);
		}
		this.config = config;
		
	}
	
	public abstract void runCmd() throws Exception;
	
	public void postRun() {}
	
	public void run(String[] args) throws ParseException, IOException{
		preRun(args);
		try {
			runCmd();
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
		postRun();
	}
	
	public String getCmdName() {
		return cmdName;
	}
	
	/**
	 * Return the options for the command.  
	 * @return Options may contain instances of ComplexOption, which specifies
	 * an alias that is used to connect Command line options with config file
	 * options.
	 */
	public Options getOptions() {
		if (this.options != null) return this.options;
		Options options = new Options();

		Option helpOption = new Option("h", "help", false, "Get Help.");
		options.addOption(helpOption);
		
		ValidatingOption configFileOption = new ValidatingOption("c", "config-file", true, "Config file");
		configFileOption.setArgs(100);
		configFileOption.addValidator(new ReadableFileValidator(configFileOption));
		configFileOption.addValidator(new MaxValuesValidator(1));
		options.addOption(configFileOption);
		
		Option printParametersOption = new Option("p", "print-parameters", false, "Print parameters and exit.");
		options.addOption(printParametersOption);
		
		addOptions(options);
		
		this.options = options;
		return this.options;
	}
	
	public abstract void addOptions(Options options);
	
}

/**
 *	Customizing org.apache.commons.cli.HelpFormatter.
 *	Fixes a problem in HelpFormatter where wrapping wasn't being done well.
 *  Uses the terminal width.
 *	Double spaces between options.
 * @author mullermw
 *
 */
class MyHelpFormatter extends HelpFormatter {
	
	public static final int TERMINAL_WIDTH;
	private static final Pattern FIRST_LINE_PATTERN;
	
	static {
		String termWidthProp = System.getProperty("term.width", "80");
		TERMINAL_WIDTH = termWidthProp.matches("^\\d+$") ? Integer.parseInt(termWidthProp) : 80;
		FIRST_LINE_PATTERN = Pattern.compile("(.{1,"+TERMINAL_WIDTH+"})(\\s|$)");
	}
	
	public MyHelpFormatter() {
		this.setWidth(TERMINAL_WIDTH);
	}
	
	@Override
	public StringBuffer renderWrappedText(StringBuffer sb, int width,
			int nextLineTabStop, String text) {

		String indent = TextUtils.getSpaces(nextLineTabStop);
		Matcher matcher = FIRST_LINE_PATTERN.matcher(text);
		if (!matcher.find()) return super.renderWrappedText(sb, width, nextLineTabStop, text);
		
		StringBuffer formatted = new StringBuffer();
		formatted.append(matcher.group(1));
		formatted.append("\n");
		text = matcher.replaceFirst("").trim();
		int descWidth = TERMINAL_WIDTH - nextLineTabStop;
		if (descWidth < 20) descWidth = 20;  //Minimum width for description is 20.
		Pattern descLinePattern = Pattern.compile("(.{1,"+descWidth+"})(\\s|$)");

		while (text.length() > descWidth) {
			matcher = descLinePattern.matcher(text);
			matcher.find();
			formatted.append(indent);
			formatted.append(matcher.group(1));
			formatted.append("\n");
			text = matcher.replaceFirst("");
		}
		if (text.length() > 0) {
			formatted.append(indent);
			formatted.append(text);
		}
		while(!formatted.toString().endsWith("\n"))
			formatted.append("\n");
		
		sb.append(formatted);

		return sb;
	}
	
	@Override
	public void printHelp(String cmdLineSyntax, Options options) {
		super.printHelp(cmdLineSyntax, options);
	}
	
	
};