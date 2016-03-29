package com.lifetechnologies.solid.wt.config.validator;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.solid.wt.mapper.Mapper;

public class EnumValidator<T extends Enum<T>> extends AbstractValidator {
	
	private Class<T> clazz;
	
	public EnumValidator(Class<T> clazz) {	
		if (clazz == null) throw new IllegalArgumentException("clazz cannot be null");
		if (clazz.isEnum() == false) throw new IllegalArgumentException(clazz.getSimpleName()+" is not an enum");
		this.clazz = clazz;
	}
	
	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		if (value == null) return new ValidatorResultImpl(false, "value is null");
		try {
			Enum.valueOf(clazz, value);
			return new ValidatorResultImpl(true);
		} catch (IllegalArgumentException e) {}
		try {
			Enum.valueOf(clazz, value.toLowerCase());
			return new ValidatorResultImpl(true);
		} catch (IllegalArgumentException e) {}
		try {
			Enum.valueOf(clazz, value.toUpperCase());
			return new ValidatorResultImpl(true);
		} catch (IllegalArgumentException e) {}
		return new ValidatorResultImpl(false, value + " is not one of " + Arrays.asList(clazz.getEnumConstants()));
		
	}

	public static void main(String[] args) {
		EnumValidator<Mapper.MappingStep> enumValidator = new EnumValidator<Mapper.MappingStep>(Mapper.MappingStep.class);
		ValidatorResult result = enumValidator.validate("foo", null);
		System.out.println(result.isValid()+ " "+result.getMessages());
		result = enumValidator.validate("splitting", null);
		System.out.println(result.isValid()+ " "+result.getMessages());
	}
}
