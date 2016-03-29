package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validates required arguments.
 * @author mullermw
 *
 */
public class RequiredArgumentValidator extends AbstractValidator {

	@Override
	public ValidatorResult validateAll(Set<String> values, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl mergedResult = new ValidatorResultImpl(true);
		if (values == null || values.size() == 0) {
			mergedResult.setValid(false);
			mergedResult.addMessage("Missing required parameter");
		} else {
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
	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {

		ValidatorResultImpl result = new ValidatorResultImpl(true);
		if (Utilities.isBlank(value)) {
			result.setValid(false);
			result.addMessage("Missing required parameter");
		}
		return result;
	}
}