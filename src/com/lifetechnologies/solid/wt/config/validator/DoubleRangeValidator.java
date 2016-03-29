package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates that the values are within a range of doubles.
 * @author mullermw
 *
 */
public class DoubleRangeValidator extends AbstractValidator {
	private Option option;
	private Double min = Double.MIN_VALUE;
	private Double max = Double.MAX_VALUE;

	public DoubleRangeValidator(Option option, double min, double max) {
		this.option = option;
		this.min = min;
		this.max = max;
	}

	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		if (value == null) value = "";
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		Double doubleValue;
		String errMsg = option.getLongOpt() + " must be a number between " + min + " & " + max + ".";
		try {
			doubleValue = Double.valueOf(value);
		} catch (NumberFormatException e) {
			result.setValid(false);
			result.addMessage(errMsg);
			return result;
		}
		if (doubleValue < min || doubleValue > max) {
			result.setValid(false);
			result.addMessage(errMsg);
			return result;
		}
		return result;
	}
}
