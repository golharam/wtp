package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validate multiple values
 * @author mullermw
 *
 */
public abstract class AbstractValidator implements Validator {
	@Override
	public ValidatorResult validateAll(Set<String> values, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl mergedResult = new ValidatorResultImpl(true);
		if (values != null ) {
			for (String value : values) {
				ValidatorResult result = this.validate(value, config);
				if (!result.isValid()) {
					mergedResult.setValid(false);
					mergedResult.addMessages(result.getMessages());
				}
			}
		}
		return mergedResult;
	}

	abstract public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config);
}
