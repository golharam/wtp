package com.lifetechnologies.solid.wt.config.validator;

import java.io.File;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

public class CreateableFileValidator extends AbstractValidator {

	private Option option;
	
	public CreateableFileValidator(Option option) {
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
		File file = new File(s).getAbsoluteFile();

		if (file.exists() && file.isDirectory()) {
			result.setValid(false);
			result.addMessage(s + " is an existing directory.");
		}
		File parentFile = file.getParentFile();
		if (!parentFile.exists()) {
			result.setValid(false);
			result.addMessage("Directory: "+parentFile+" does not exist.");
			return result;
		}
		if (!parentFile.canWrite()) {
			result.setValid(false);
			result.addMessage("File: "+s+" is not writeable.");
		}
		if (file.exists() && !file.canWrite()) {
			result.setValid(false);
			result.addMessage("File: "+s+" is not writeable.");
		}
		return result;
	}
}
