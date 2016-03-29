package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validate a set of values.
 * @author mullermw
 *
 */
public interface Validator {
	/**
	 * Validate all the values
	 * @param value to validate
	 * @return the result.
	 */
	public ValidatorResult validateAll(Set<String> value, Map<ConfigKey, Set<String>> config);
}
