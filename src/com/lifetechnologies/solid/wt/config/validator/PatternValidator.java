package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates that the values match a pattern.
 * @author mullermw
 *
 */
public class PatternValidator extends AbstractValidator {

	Option option;
	Pattern pattern;

	public PatternValidator(Option option, Pattern pattern) {
		this.option = option;
		this.pattern = pattern;
	}

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		if (pattern == null) return result;
		if (value == null) value = "";
		if (!pattern.matcher(value).matches()) {
			result.setValid(false);
			result.addMessage("Invalid value for "+option.getLongOpt()+": '"+value+"'");
		}
		return result;
	}
}
