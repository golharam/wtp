package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates the number of arguments.
 * @author mullermw
 *
 */
public class MaxValuesValidator implements Validator {
	
	int maxValues = 0;
	
	public MaxValuesValidator(int maxValues) {
		this.maxValues = maxValues;
	}
	
	@Override
	public ValidatorResult validateAll(Set<String> values, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		if (values != null && values.size() > maxValues ) {
			result.setValid(false);
			result.addMessage("Too many values.");
		}
		return result;
	}
}
