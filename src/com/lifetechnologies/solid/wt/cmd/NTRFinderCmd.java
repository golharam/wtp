package com.lifetechnologies.solid.wt.cmd;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.regex.Pattern;

import org.apache.commons.cli.Options;
import com.lifetechnologies.solid.wt.NTRFinderLauncher;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.ComplexOption;
import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.solid.wt.config.validator.CreatableDirectoryValidator;
import com.lifetechnologies.solid.wt.config.validator.DoubleRangeValidator;
import com.lifetechnologies.solid.wt.config.validator.GenomicRegionValidator;
import com.lifetechnologies.solid.wt.config.validator.IntegerRangeValidator;
import com.lifetechnologies.solid.wt.config.validator.MaxValuesValidator;
import com.lifetechnologies.solid.wt.config.validator.PatternValidator;
import com.lifetechnologies.solid.wt.config.validator.ReadableFileValidator;
import com.lifetechnologies.solid.wt.config.validator.RequiredArgumentValidator;

public class NTRFinderCmd extends AbstractCmd {

	private static final String NAME = "NTR Finder";
	
	public NTRFinderCmd() {
		super(NAME);
	}
	
	/**
	 * Main entry point for the NTR Finder.
	 * @param args are parsed as command line switches.
	 */
	public void runCmd() {
		try {
			
			File outputFolder = new File(config.get(ConfigKey.WT_OUTPUT_DIR).iterator().next());
			if (!outputFolder.exists()) {
				outputFolder.mkdir();
				if (!outputFolder.exists()) throw new Exception("Cannot open output dir: "+outputFolder);
			}
			if (outputFolder.isDirectory() == false) throw new Exception(outputFolder +" is a directory.");
			if (outputFolder.canWrite() == false) throw new Exception("Cannot write to output folder: "+outputFolder);

			
			File logFile = new File(outputFolder, "ntr_finder.log");
			if (logFile.exists()) {
				//Rename old log file.
				int logFileCounter = 1;
				File oldLogFile = null;
				do {
					oldLogFile = new File(outputFolder, "ntr_finder."+(logFileCounter++)+".log");
				} while (oldLogFile.exists());
				logFile.renameTo(oldLogFile);
			}

			Handler handler = new FileHandler(logFile.getPath());
			handler.setFormatter(new SimpleFormatter());
			Logger rootLogger = LogManager.getLogManager().getLogger("");
			rootLogger.addHandler(handler);
			System.err.println("Logging messages to "+logFile);
			
			//Log config
			List<ConfigKey> configKeys = new ArrayList<ConfigKey>(config.keySet());
			Collections.sort(configKeys);
			StringBuffer configBuffer = new StringBuffer("Configuration values: \n");
			for (ConfigKey key : configKeys ) {
				List<String> values = new ArrayList<String>(config.get(key));
				List<Integer> numValues = null;
				try {
					numValues = Utilities.toIntegerList(values);
					Collections.sort(numValues);
				} catch (NumberFormatException e) {}
				Collections.sort(values);
				if (numValues != null) {
					for (Integer value : numValues)
						configBuffer.append(key.getQualifiedName().concat("=").concat(value.toString()).concat("\n"));
				} else {
					for (String value : values)
						configBuffer.append(key.getQualifiedName().concat("=").concat(value).concat("\n"));
				}
			}
			logger.info(configBuffer.toString());
			
			Set<String> filenameAtrReference = config.get(ConfigKey.WT_NTR_FILE_ATR_REFERENCE);
			if (filenameAtrReference == null || filenameAtrReference.isEmpty())
				logger.warning("No exon reference file specified (file-atr-reference-gff).");

			//Config OK, proceed.
			boolean submitJobs = true;
			if (config.containsKey(ConfigKey.WT_SKIP_JOBS)) {
				if (config.get(ConfigKey.WT_SKIP_JOBS) == null ||
					config.get(ConfigKey.WT_SKIP_JOBS).size() == 0) {
					submitJobs = false;
				} else if (Boolean.parseBoolean(config.get(ConfigKey.WT_SKIP_JOBS).iterator().next())) {
					submitJobs = false;
				}
			}
			logger.fine("submitJobs="+submitJobs);
			new NTRFinderLauncher().launch(config, submitJobs);
			
			logger.info("NTR Finding completed successfully.");
		} catch (Exception e) {
			logger.severe("An error has stopped NTR Finding.\n" + Utilities.toStackTrace(e));
			System.exit(1);
		}
		System.err.println("Finished.");
		logger.info("Exit");
		System.exit(0);
	}

	public void addOptions(Options options) {

		// Note that for all options that accept arguments, I setArgs(100).  This
		// is because I am checking the number of arguments with MaxValuesValidator.
		ComplexOption schedulingEnvironmentOption = new ComplexOption(ConfigKey.QUEUE_SYS, true);
		schedulingEnvironmentOption.setArgs(100);
		schedulingEnvironmentOption.addValidator(new RequiredArgumentValidator());
		schedulingEnvironmentOption.addValidator(new MaxValuesValidator(1));
		schedulingEnvironmentOption.addValidator(new PatternValidator(schedulingEnvironmentOption, Pattern.compile("^pbs|lsf|sge$", Pattern.CASE_INSENSITIVE)));
		options.addOption(schedulingEnvironmentOption);
		
		ComplexOption queueNameOption = new ComplexOption(ConfigKey.QUEUE_SYS_QUEUE, true);
		queueNameOption.setArgs(100);
		queueNameOption.addValidator(new RequiredArgumentValidator());
		queueNameOption.addValidator(new MaxValuesValidator(1));
		queueNameOption.addValidator(new PatternValidator(queueNameOption, Pattern.compile("\\S+")));
		options.addOption(queueNameOption);
		
		ComplexOption schedResReqOption = new ComplexOption(ConfigKey.QUEUE_SYS_RESOURCE_STRING, true);
		schedResReqOption.setArgs(100);
		schedResReqOption.addValidator(new MaxValuesValidator(1));
		options.addOption(schedResReqOption);
		
		ComplexOption additionSchedOption = new ComplexOption(ConfigKey.QUEUE_SYS_OPTIONS, true);
		additionSchedOption.setArgs(100);
		additionSchedOption.addValidator(new MaxValuesValidator(1));
		options.addOption(additionSchedOption);
		
		ComplexOption maxFileOption = new ComplexOption(ConfigKey.WT_NTR_MAX_FILE, true);
		maxFileOption.setArgs(100);
		maxFileOption.setValueSeparator(',');
		maxFileOption.addValidator(new RequiredArgumentValidator());
		maxFileOption.addValidator(new ReadableFileValidator(maxFileOption));
		options.addOption(maxFileOption);

		ComplexOption windowSizeOption = new ComplexOption(ConfigKey.WT_NTR_MIN_WINDOW_SIZE,  true);
		windowSizeOption.setArgs(100);
		windowSizeOption.setValueSeparator(',');
		windowSizeOption.addValidator(new RequiredArgumentValidator());
		windowSizeOption.addValidator(new IntegerRangeValidator(windowSizeOption, 25, 10000));
		options.addOption(windowSizeOption);

		ComplexOption minWindowCoverageOption = new ComplexOption(ConfigKey.WT_NTR_MIN_WINDOW_COVERAGE, true);
		minWindowCoverageOption.setArgs(100);
		minWindowCoverageOption.setValueSeparator(',');
		minWindowCoverageOption.addValidator(new RequiredArgumentValidator());
		minWindowCoverageOption.addValidator(new DoubleRangeValidator(minWindowCoverageOption, 0.01, 100000));
		options.addOption(minWindowCoverageOption);

		ComplexOption overlapOption = new ComplexOption(ConfigKey.WT_NTR_MIN_OVERLAP, true);
		overlapOption.setArgs(100);
		overlapOption.addValidator(new RequiredArgumentValidator());
		overlapOption.addValidator(new MaxValuesValidator(1));
		overlapOption.addValidator(new DoubleRangeValidator(overlapOption, 0.001, 1.0));
		options.addOption(overlapOption);

		ComplexOption minScoreOption = new ComplexOption(ConfigKey.WT_NTR_MIN_ALIGNMENT_SCORE, true);
		minScoreOption.setArgs(100);
		minScoreOption.addValidator(new RequiredArgumentValidator());
		minScoreOption.addValidator(new MaxValuesValidator(1));
		minScoreOption.addValidator(new IntegerRangeValidator(minScoreOption, 1, Integer.MAX_VALUE));
		options.addOption(minScoreOption);

		ComplexOption trimmingFractionOption = new ComplexOption(ConfigKey.WT_NTR_TRIMMING_FRACTION, true);
		trimmingFractionOption.setArgs(100);
		trimmingFractionOption.addValidator(new RequiredArgumentValidator());
		trimmingFractionOption.addValidator(new MaxValuesValidator(1));
		trimmingFractionOption.addValidator(new DoubleRangeValidator(trimmingFractionOption, 0.0, 1.0));
		options.addOption(trimmingFractionOption);

		ComplexOption genomicRegionsOption = new ComplexOption(ConfigKey.WT_NTR_GENOMIC_REGION, true);
		genomicRegionsOption.setArgs(100);
		genomicRegionsOption.setValueSeparator(',');
		genomicRegionsOption.addValidator(new GenomicRegionValidator(genomicRegionsOption));
		options.addOption(genomicRegionsOption);
		
		ComplexOption maxNTROption = new ComplexOption(ConfigKey.WT_NTR_MAX_PTRS_PER_MEGABASE, true);
		maxNTROption.setArgs(100);
		maxNTROption.addValidator(new RequiredArgumentValidator());
		maxNTROption.addValidator(new MaxValuesValidator(1));
		maxNTROption.addValidator(new IntegerRangeValidator(maxNTROption, 0, Integer.MAX_VALUE));
		options.addOption(maxNTROption);
		
		ComplexOption fileReferenceOption = new ComplexOption(ConfigKey.WT_FILE_REFERENCE, true);
		fileReferenceOption.setArgs(100);
		fileReferenceOption.addValidator(new RequiredArgumentValidator());
		fileReferenceOption.addValidator(new MaxValuesValidator(1));
		fileReferenceOption.addValidator(new ReadableFileValidator(fileReferenceOption));
		options.addOption(fileReferenceOption);
		
		ComplexOption fileExonReferenceOption = new ComplexOption(ConfigKey.WT_NTR_FILE_ATR_REFERENCE, true );
		fileExonReferenceOption.setArgs(100);
		fileExonReferenceOption.addValidator(new MaxValuesValidator(1));
		fileExonReferenceOption.addValidator(new ReadableFileValidator(fileExonReferenceOption));
		options.addOption(fileExonReferenceOption);

		ComplexOption outputDirOption = new ComplexOption(ConfigKey.WT_OUTPUT_DIR, true );
		outputDirOption.setArgs(100);
		outputDirOption.addValidator(new RequiredArgumentValidator());
		outputDirOption.addValidator(new MaxValuesValidator(1));
		outputDirOption.addValidator(new CreatableDirectoryValidator(outputDirOption));
		options.addOption(outputDirOption);
		
		ComplexOption skipJobSubOption = new ComplexOption(ConfigKey.WT_SKIP_JOBS, false );
		options.addOption(skipJobSubOption);
		
		ComplexOption deleteIntermediateFilesOption = new ComplexOption(ConfigKey.WT_DELETE_INTERMEDIATE_FILES, false);
		options.addOption(deleteIntermediateFilesOption);
		
	}
	
}





























