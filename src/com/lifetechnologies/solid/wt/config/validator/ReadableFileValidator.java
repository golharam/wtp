package com.lifetechnologies.solid.wt.config.validator;

import java.io.File;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates that the values are readable files.
 * @author mullermw
 *
 */
public class ReadableFileValidator extends AbstractValidator {

	Option option;

	public ReadableFileValidator(Option option) {
		this.option = option;
	}

	@Override
	public ValidatorResult validate(String s, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		if (s == null) {
			result.setValid(false);
			result.addMessage(option.getLongOpt() + " Cannot be null.");
			return result;
		}
		if (s.isEmpty()) {
			result.setValid(false);
			result.addMessage(option.getLongOpt() + " Cannot be empty.");
			return result;
		}
		File file = new File(s);
		if (!file.exists()) {
			result.setValid(false);
			result.addMessage("File: "+s+" does not exist.");
			return result;
		}
		if (!file.isFile()) {
			result.setValid(false);
			result.addMessage("File: "+s+" is not a file.");
		}
		if (!file.canRead()) {
			result.setValid(false);
			result.addMessage("File: "+s+" is not readable.");
		}
		return result;
	}
}
