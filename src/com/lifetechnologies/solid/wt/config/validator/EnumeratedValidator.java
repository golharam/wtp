package com.lifetechnologies.solid.wt.config.validator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;


public class EnumeratedValidator extends AbstractValidator {

	Set<String> values = new HashSet<String>();
	
 	public EnumeratedValidator(String ...strings ) {
 		values.addAll(Arrays.asList(strings));
	}
	
	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		if (value == null) return new ValidatorResultImpl(false, "value is null.");
		if (values.contains(value.toLowerCase().trim()))return new ValidatorResultImpl(true);
		return new ValidatorResultImpl(false, "'"+value+"' is not one of the legal values: "+values);
	}

}
