package com.lifetechnologies.solid.wt.config;

import com.lifetechnologies.solid.wt.config.validator.MaxValuesValidator;
import com.lifetechnologies.solid.wt.config.validator.RequiredArgumentValidator;


/**
 * Extends Option with an third Option identifier, alias.  I'm using this
 * to connect the command line options with config file options.
 * @author mullermw
 *
 */
public class ComplexOption extends ValidatingOption {

	public static final long serialVersionUID = 1;
	
	private ConfigKey key;
	
	public ComplexOption(ConfigKey key, boolean hasArg) 
		throws IllegalArgumentException {
		super(key.asCommandLineShortOpt(), key.asCommandLineLongOpt(), hasArg, key.getDescription());
		this.key = key;
	}

	public ConfigKey getKey() {
		return key;
	}
	
	public static ComplexOption newRequiredSingularOption(ConfigKey key, boolean hasArg) {
		ComplexOption option = new ComplexOption(key, hasArg);
		option.setArgs(100);
		option.addValidator(new RequiredArgumentValidator());
		option.addValidator(new MaxValuesValidator(1));
		return option;
	}
	
	public static ComplexOption newOptionalSingularOption(ConfigKey key, boolean hasArg) {
		ComplexOption option = new ComplexOption(key, hasArg);
		option.setArgs(100);
		option.addValidator(new MaxValuesValidator(1));
		return option;
	}
	
}
