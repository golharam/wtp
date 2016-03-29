package com.lifetechnologies.solid.wt.cmd;

import java.io.File;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.regex.Pattern;

import org.apache.commons.cli.Options;

import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.ComplexOption;
import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.solid.wt.config.validator.BooleanValidator;
import com.lifetechnologies.solid.wt.config.validator.ByteValueValidator;
import com.lifetechnologies.solid.wt.config.validator.CreatableDirectoryValidator;
import com.lifetechnologies.solid.wt.config.validator.DoubleRangeValidator;
import com.lifetechnologies.solid.wt.config.validator.EnumValidator;
import com.lifetechnologies.solid.wt.config.validator.EnumeratedValidator;
import com.lifetechnologies.solid.wt.config.validator.IntegerRangeValidator;
import com.lifetechnologies.solid.wt.config.validator.PatternValidator;
import com.lifetechnologies.solid.wt.config.validator.ReadMaskValidator;
import com.lifetechnologies.solid.wt.config.validator.ReadableFileValidator;
import com.lifetechnologies.solid.wt.config.validator.Validator;
import com.lifetechnologies.solid.wt.config.validator.ValidatorResult;
import com.lifetechnologies.solid.wt.config.validator.ValidatorResultImpl;
import com.lifetechnologies.solid.wt.mapper.FilteringMode;
import com.lifetechnologies.solid.wt.mapper.Mapper;
import com.lifetechnologies.util.MathUtils;

import static com.lifetechnologies.solid.wt.Utilities.firstValue;

/**
 * Maps reads to reference sequence.
 * @author mullermw
 *
 */
public class MapperCmd extends AbstractCmd {

	private static final String NAME = "Split Read Mapper";
	
	public MapperCmd() {
		super(NAME);
	}
	
	@Override
	public void addOptions(Options options) {

		ComplexOption schedulingEnvironmentOption = ComplexOption.newRequiredSingularOption(ConfigKey.QUEUE_SYS, true);
		schedulingEnvironmentOption.addValidator(new EnumeratedValidator("pbs", "lsf", "sge"));
		options.addOption(schedulingEnvironmentOption);
		
		ComplexOption queueNameOption = ComplexOption.newRequiredSingularOption(ConfigKey.QUEUE_SYS_QUEUE, true);
		queueNameOption.addValidator(new PatternValidator(queueNameOption, Pattern.compile("\\S+")));
		options.addOption(queueNameOption);
		
		ComplexOption schedResReqOption = ComplexOption.newOptionalSingularOption(ConfigKey.QUEUE_SYS_RESOURCE_STRING, true);
		options.addOption(schedResReqOption);
		
		ComplexOption additionSchedOption = ComplexOption.newOptionalSingularOption(ConfigKey.QUEUE_SYS_OPTIONS, true);
		options.addOption(additionSchedOption);
		
		ComplexOption tmpDirOption = ComplexOption.newRequiredSingularOption(ConfigKey.QUEUE_SYS_NODES_TMP_DIR, true);
		options.addOption(tmpDirOption);
		
		ComplexOption maxMemoryPerJobOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAX_MEMORY_PER_JOB, true);
		maxMemoryPerJobOption.addValidator(new ByteValueValidator());
		options.addOption(maxMemoryPerJobOption);
		
		ComplexOption memoryRequirementAdjustmentOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR, true);
		memoryRequirementAdjustmentOption.addValidator(new DoubleRangeValidator(memoryRequirementAdjustmentOption, 0.0, Double.MAX_VALUE));
		options.addOption(memoryRequirementAdjustmentOption);
		
		ComplexOption deleteIntermediateFilesOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_DELETE_INTERMEDIATE_FILES, true);
		deleteIntermediateFilesOption.addValidator(new BooleanValidator());
		options.addOption(deleteIntermediateFilesOption);
		
		ComplexOption fileReferenceOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_FILE_REFERENCE, true);
		fileReferenceOption.addValidator(new ReadableFileValidator(fileReferenceOption));
		options.addOption(fileReferenceOption);
		
		ComplexOption fileExonReferenceOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_EXON_REFERENCE, true);
		fileExonReferenceOption.addValidator(new ReadableFileValidator(fileExonReferenceOption));
		options.addOption(fileExonReferenceOption);
		
		ComplexOption outputDirOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_OUTPUT_DIR, true);
		outputDirOption.addValidator(new CreatableDirectoryValidator(outputDirOption));
		options.addOption(outputDirOption);
		
		ComplexOption skipToOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_TASK_SKIP_TO, true);
		skipToOption.addValidator(new EnumValidator<Mapper.MappingStep>(Mapper.MappingStep.class));
		options.addOption(skipToOption);
		
		ComplexOption fileFilterReferenceOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_FILTER_FILE_REFERENCE, true);
		fileFilterReferenceOption.addValidator(new ReadableFileValidator(fileFilterReferenceOption));
		options.addOption(fileFilterReferenceOption);
		
		ComplexOption filterModeOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_FILTER_MODE, true);
		filterModeOption.addValidator(new EnumValidator<FilteringMode>(FilteringMode.class));
		options.addOption(filterModeOption);
		
		ComplexOption filterMismatchOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_FILTER_MISMATCHES, true);
		filterMismatchOption.addValidator(new Validator() {
			@Override
			public ValidatorResult validateAll(Set<String> value,
					Map<ConfigKey, Set<String>> config) {
				String filterMode = firstValue(config.get(ConfigKey.WT_MAPPING_FILTER_MODE));
				if (filterMode == null || filterMode.toLowerCase().trim().equals("off")) return new ValidatorResultImpl(true);
				if ((value == null || value.isEmpty()) &&
					(filterMode.toLowerCase().trim().equals("one_or_more") ||
					 filterMode.toLowerCase().trim().equals("both"))) {
					return new ValidatorResultImpl(false, ConfigKey.WT_MAPPING_FILTER_MISMATCHES.getQualifiedName() + " is required when filtering is on.");
				}
				return new ValidatorResultImpl(true);
			}
		});
		filterMismatchOption.addValidator(new IntegerRangeValidator(filterMismatchOption, 0, 10));
		options.addOption(filterMismatchOption);
				
		ComplexOption readMaskOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_READ_MASK, true);
		readMaskOption.addValidator(new ReadMaskValidator());
		options.addOption(readMaskOption);
		
		ComplexOption maxLociOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_MAX_LOCI, true);
		maxLociOption.addValidator(new IntegerRangeValidator(maxLociOption, 1,100));
		options.addOption(maxLociOption);
		
		ComplexOption validAdjacentOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE, true);
		validAdjacentOption.addValidator(new BooleanValidator());
		options.addOption(validAdjacentOption);
		
		ComplexOption iubOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_MATCH_IUB, true);
		iubOption.addValidator(new BooleanValidator());
		options.addOption(iubOption);
		
		ComplexOption readLengthOption = ComplexOption.newOptionalSingularOption(ConfigKey.WT_MAPPING_READ_LENGTH, true);
		readLengthOption.addValidator(new IntegerRangeValidator(readLengthOption, 25, 200));
		options.addOption(readLengthOption);
		
		ComplexOption leftLengthOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SPLIT_LEFT_LENGTH, true);
		leftLengthOption.addValidator(new IntegerRangeValidator(leftLengthOption, 0, 100));
		leftLengthOption.addValidator(new NotGreaterThanReadLengthValidator());
		options.addOption(leftLengthOption);
		
		ComplexOption leftMismatchesOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SPLIT_LEFT_MISMATCHES, true);
		leftMismatchesOption.addValidator(new IntegerRangeValidator(leftMismatchesOption, 0,100));
		leftMismatchesOption.addValidator(new NotGreaterThanReadLengthValidator());
		options.addOption(leftMismatchesOption);
		
		ComplexOption rightLengthOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SPLIT_RIGHT_LENGTH, true);
		rightLengthOption.addValidator(new IntegerRangeValidator(rightLengthOption, 0, 100));
		rightLengthOption.addValidator(new NotGreaterThanReadLengthValidator());
		options.addOption(rightLengthOption);
		
		ComplexOption rightMismatchesOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SPLIT_RIGHT_MISMATCHES, true);
		rightMismatchesOption.addValidator(new IntegerRangeValidator(rightMismatchesOption, 0, 100));
		rightMismatchesOption.addValidator(new NotGreaterThanReadLengthValidator());
		options.addOption(rightMismatchesOption);
		
		ComplexOption minScoreOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SCORE_MIN, true);
		minScoreOption.addValidator(new IntegerRangeValidator(minScoreOption, 0, 100));
		options.addOption(minScoreOption);
		
		ComplexOption uniqenessGapOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_SCORE_UNIQUENESS_GAP, true);
		uniqenessGapOption.addValidator(new IntegerRangeValidator(uniqenessGapOption, 0, 100));
		options.addOption(uniqenessGapOption);
		
		ComplexOption readsFileOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_FILE_READS, true);
		readsFileOption.addValidator(new ReadableFileValidator(readsFileOption));
		options.addOption(readsFileOption);
		
		ComplexOption penaltyKnownJunctionOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_PENALTY_KNOWN_JUNCTION, true);
		penaltyKnownJunctionOption.addValidator(new IntegerRangeValidator(penaltyKnownJunctionOption, 0, 50));
		options.addOption(penaltyKnownJunctionOption);
		
		ComplexOption penaltyPutativeJunctionOption = ComplexOption.newRequiredSingularOption(ConfigKey.WT_MAPPING_PENALTY_PUTATIVE_JUNCTION, true);
		penaltyPutativeJunctionOption.addValidator(new IntegerRangeValidator(penaltyPutativeJunctionOption, 0, 50));
		options.addOption(penaltyPutativeJunctionOption);
		
	}
	
	@Override
	public void runCmd() throws Exception {
		
		File outputFolder = new File(config.get(ConfigKey.WT_OUTPUT_DIR).iterator().next());
		if (!outputFolder.exists()) {
			outputFolder.mkdir();
			if (!outputFolder.exists()) throw new Exception("Cannot open output dir: "+outputFolder);
		}
		if (outputFolder.isDirectory() == false) throw new Exception(outputFolder +" is a directory.");
		if (outputFolder.canWrite() == false) throw new Exception("Cannot write to output folder: "+outputFolder);

		
		File logFile = new File(outputFolder, "mapper.log");
		if (logFile.exists()) {
			//Rename old log file.
			int logFileCounter = 1;
			File oldLogFile = null;
			do {
				oldLogFile = new File(outputFolder, "mapper."+(logFileCounter++)+".log");
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
		
		Mapper mapper = new Mapper();
		writeConfigTo(mapper);
		mapper.mapReads();
	}
	
	private void writeConfigTo(Mapper mapper) throws ParseException {
		
		mapper.setSchedulingEnvironment(first(ConfigKey.QUEUE_SYS));
		mapper.setNameOfQueue(first(ConfigKey.QUEUE_SYS_QUEUE));
		mapper.setMaxMemoryPerJob(MathUtils.toBytes(first(ConfigKey.WT_MAX_MEMORY_PER_JOB)));
		mapper.setMemoryRequirementAdjustmentFactor(firstDouble(ConfigKey.WT_MAPPING_MEMORY_REQUIREMENT_ADJUSTMENT_FACTOR));
		mapper.setSchedulerResourceRequirements(first(ConfigKey.QUEUE_SYS_RESOURCE_STRING));
		mapper.setAdditionalSchedulerOptions(first(ConfigKey.QUEUE_SYS_OPTIONS));
		mapper.setSkipTo(Mapper.MappingStep.parseMappingStep(first(ConfigKey.WT_MAPPING_TASK_SKIP_TO)));
		mapper.setFilteringMode(FilteringMode.parseFilteringMode(first(ConfigKey.WT_MAPPING_FILTER_MODE)));
		mapper.setReadLength(firstInt(ConfigKey.WT_MAPPING_READ_LENGTH));
		mapper.setReadMask(new ReadMask(first(ConfigKey.WT_MAPPING_READ_MASK)));
		mapper.setMaxMappingLocations(firstInt(ConfigKey.WT_MAPPING_MAX_LOCI));
		mapper.setValidAdjacentRules(firstBoolean(ConfigKey.WT_MAPPING_VALID_ADJACENT_MISMATCHES_COUNT_AS_ONE));
		mapper.setMatchIub(firstBoolean(ConfigKey.WT_MAPPING_MATCH_IUB));
		mapper.setLeftLength(firstInt(ConfigKey.WT_MAPPING_SPLIT_LEFT_LENGTH));
		mapper.setRightLength(firstInt(ConfigKey.WT_MAPPING_SPLIT_RIGHT_LENGTH));
		mapper.setLeftMismatches(firstInt(ConfigKey.WT_MAPPING_SPLIT_LEFT_MISMATCHES));
		mapper.setRightMismatches(firstInt(ConfigKey.WT_MAPPING_SPLIT_RIGHT_MISMATCHES));
		mapper.setFilterMismatches(firstInt(ConfigKey.WT_MAPPING_FILTER_MISMATCHES));
		mapper.setMinScore(firstInt(ConfigKey.WT_MAPPING_SCORE_MIN));
		mapper.setScoreGap(firstInt(ConfigKey.WT_MAPPING_SCORE_UNIQUENESS_GAP));
		mapper.setFilterReference(firstFile(ConfigKey.WT_MAPPING_FILTER_FILE_REFERENCE));
		mapper.setReferenceFile(firstFile(ConfigKey.WT_FILE_REFERENCE));
		mapper.setExonReference(firstFile(ConfigKey.WT_EXON_REFERENCE));
		mapper.setCsfastaFile(firstFile(ConfigKey.WT_MAPPING_FILE_READS));
		mapper.setScratchDir(firstFile(ConfigKey.QUEUE_SYS_NODES_TMP_DIR));
		mapper.setOutputDir(firstFile(ConfigKey.WT_OUTPUT_DIR));
		mapper.setKnownJunctionPenalty(firstInt(ConfigKey.WT_MAPPING_PENALTY_KNOWN_JUNCTION));
		mapper.setPutativeJunctionPenalty(firstInt(ConfigKey.WT_MAPPING_PENALTY_PUTATIVE_JUNCTION));
		mapper.setDeleteIntermediateFiles(firstBoolean(ConfigKey.WT_DELETE_INTERMEDIATE_FILES));
	}
	
	private String first(ConfigKey key) {
		String val = firstValue(config.get(key));
		return val;
	}
	
	private int firstInt(ConfigKey key) {
		return Integer.parseInt(first(key));
	}
	
	private double firstDouble(ConfigKey key) {
		return Double.parseDouble(first(key));
	}
	
	private boolean firstBoolean(ConfigKey key) {
		return Utilities.isTrue(first(key));
	}
	
	private File firstFile(ConfigKey key) {
		String val = first(key);
		if (val == null) return null;
		return new File(val);
	}
	
	class NotGreaterThanReadLengthValidator implements Validator {
		@Override
		public ValidatorResult validateAll(Set<String> values,
				Map<ConfigKey, Set<String>> config) {
			if (values == null) values = new HashSet<String>();
			try {
				for (String value : values) {
					int readLength = Integer.parseInt(firstValue(config.get(ConfigKey.WT_MAPPING_READ_LENGTH)));
					int leftLength = Integer.parseInt(value);
					if (leftLength > readLength) throw new Exception();
				}
			} catch (Exception e) {
				return new ValidatorResultImpl(false, "must be less than "+ConfigKey.WT_MAPPING_READ_LENGTH.getQualifiedName());
			}
			return new ValidatorResultImpl(true);
		}
	};
}
