package com.lifetechnologies.solid.wt.config;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.config.validator.Validator;
import com.lifetechnologies.solid.wt.config.validator.ValidatorResult;
import com.lifetechnologies.solid.wt.config.validator.ValidatorResultImpl;

/**
 * Option with a Validator used to validate the option value.
 * @author mullermw
 *
 */
public class ValidatingOption extends Option implements Validator {
	
	public static final long serialVersionUID = 1;
	
	private Set<Validator> validators = new HashSet<Validator>();

	public ValidatingOption(String opt, boolean hasArg, String description)
	throws IllegalArgumentException {
		super(opt, hasArg, description);
	}

	public ValidatingOption(String opt, String longOpt, boolean hasArg,
			String description) throws IllegalArgumentException {
		super(opt, longOpt, hasArg, description);
	}

	public ValidatingOption(String opt, String description)
	throws IllegalArgumentException {
		super(opt, description);
	}

	public void addValidator(Validator validator) {
		if (validator == null) return;
		this.validators.add(validator);
	}

	public ValidatorResult validateAll(String[] values, Map<ConfigKey, Set<String>> config) {
		Set<String> set = new HashSet<String>();
		set.addAll(Arrays.asList(values));
		return validateAll(set, config);
	}
	
	@Override
	public ValidatorResult validateAll(Set<String> values, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl masterResult = new ValidatorResultImpl(true);
		for (Validator validator : validators) {
			ValidatorResult result = validator.validateAll(values, config);
			if (!result.isValid()) {
				masterResult.setValid(false);
				masterResult.addMessages(result.getMessages());
			}
		}
		return masterResult;
	}

}
