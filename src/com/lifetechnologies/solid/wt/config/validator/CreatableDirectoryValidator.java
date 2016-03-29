package com.lifetechnologies.solid.wt.config.validator;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates that the values are Createable directories.
 * @author mullermw
 *
 */
public class CreatableDirectoryValidator extends AbstractValidator {

	Option option;

	public CreatableDirectoryValidator(Option option) {
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

		if (file.exists() && !file.isDirectory()) {
			result.setValid(false);
			result.addMessage(s + " is not a directory.");
		}
		File parentFile = file.getParentFile();
		if (!parentFile.exists()) {
			result.setValid(false);
			result.addMessage("Directory: "+parentFile+" does not exist.");
			return result;
		}
		if (!parentFile.canWrite()) {
			result.setValid(false);
			result.addMessage("Directory: "+s+" is not writeable.");
		}
		return result;
	}
	public static void main(String[] args) throws IOException {
		System.out.println(new File("ignore").getAbsoluteFile().getParent());
	}
	
}
