package com.lifetechnologies.solid.wt.config.validator;

import java.text.ParseException;
import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.util.MathUtils;

/**
 * Validate that a string can be parsed as a quantity of bytes.
 * @author mullermw
 *
 */
public class ByteValueValidator extends AbstractValidator {

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		ValidatorResult result = null;
		try {
			MathUtils.toBytes(value);
			result = new ValidatorResultImpl(true);
		} catch (ParseException e) {
			result = new ValidatorResultImpl(false, value + " Cannot be parsed as a Byte String.");
		}
		return result;
	}
	
	
}
