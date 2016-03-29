package com.lifetechnologies.solid.wt.config.validator;

import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Validate that the values can be parsed as Genomic Regions.
 * @author mullermw
 *
 */
public class GenomicRegionValidator extends AbstractValidator {

	private Option option;

	public GenomicRegionValidator(Option option) {
		this.option = option;
	}

	@Override
	public ValidatorResult validate(String value, Map<ConfigKey, Set<String>> config) {
		ValidatorResultImpl result = new ValidatorResultImpl(true);
		if (Utilities.isBlank(value)) return result;
		Set<ContiguousGenomicFeature> features = ContiguousGenomicFeature.parseFeatures(value);
		if (features.size() < 1 ) {
			result.setValid(false);
			result.addMessage("Value of "+option.getLongOpt()+": " + value+ " does not define genomic regions.");
		}
		return result;
	}
}
