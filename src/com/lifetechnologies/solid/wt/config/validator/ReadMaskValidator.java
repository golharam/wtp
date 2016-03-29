package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.ReadMask;
import com.lifetechnologies.solid.wt.config.ConfigKey;

public class ReadMaskValidator extends AbstractValidator {

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		try {
			new ReadMask(value);
		} catch (Throwable t) {
			return new ValidatorResultImpl(false, t.getMessage());
		}
		return new ValidatorResultImpl(true);
	}

}
