package com.lifetechnologies.solid.wt.config.validator;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates that the values are within a range of integers.
 * @author mullermw
 *
 */
public class IntegerRangeValidator extends AbstractValidator{
	private Option option;
	private Integer min = Integer.MIN_VALUE;
	private Integer max = Integer.MAX_VALUE;

	public IntegerRangeValidator(Option option, int min, int max) {
		this.option = option;
		this.min = min;
		this.max = max;
	}

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		final NumberFormat format = DecimalFormat.getInstance();
		if (value == null) value = "";
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		Integer intValue;
		String errMsg = "Invalid value: '"+value+"'. "+option.getLongOpt() + " must be an integer between " + min + " & " + max + ".";
		try {
			intValue = format.parse(value).intValue();
		} catch (Exception e) {
			result.setValid(false);
			result.addMessage(errMsg);
			return result;
		}
		if (intValue < min || intValue > max) {
			result.setValid(false);
			result.addMessage(errMsg);
			return result;
		}
		return result;
	}
}
