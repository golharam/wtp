package com.lifetechnologies.solid.wt.config.validator;

import java.util.List;

/**
 * The result of a Validation.
 * @author mullermw
 *
 */
public interface ValidatorResult {

	/**
	 * true if value is valid.
	 * @return
	 */
	public boolean isValid();

	/**
	 * @return a list of messages from the Validator.
	 */
	public List<String> getMessages();
}
