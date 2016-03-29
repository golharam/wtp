package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.ConfigKey;

public class BooleanValidator extends AbstractValidator {

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		if (Utilities.isTrue(value) || Utilities.isFalse(value)) {
			return new ValidatorResultImpl(true);
		} else {
			return new ValidatorResultImpl(false, "'"+value+"' cannot be parsed as a boolean." );
		}
		
	}
	
}
