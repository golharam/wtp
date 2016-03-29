package com.lifetechnologies.solid.wt.config.validator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Simple implementation of ValidatorResult.
 * @author mullermw
 *
 */
public class ValidatorResultImpl implements ValidatorResult {

	private boolean valid;
	private ArrayList<String> messages = new ArrayList<String>();

	public ValidatorResultImpl() {
	}

	public ValidatorResultImpl(boolean valid) {
		this.valid = valid;
	}

	public ValidatorResultImpl(boolean valid, String message) {
		this.valid = valid;
		this.addMessage(message);
	}

	public ValidatorResultImpl(boolean valid, Collection<String> messages) {
		this.valid = valid;
		this.addMessages(messages);
	}

	@Override
	public boolean isValid() {
		return valid;
	}

	public void setValid(boolean valid) {
		this.valid = valid;
	}

	@SuppressWarnings("unchecked")
	public List<String> getMessages() {
		return (List<String>)messages.clone();
	}

	public void addMessage(String message) {
		if (message == null) return;
		this.messages.add(message);
	}

	public void addMessages(Collection<String> messages) {
		if (messages == null) return;
		this.messages.addAll(messages);
	}	
}
