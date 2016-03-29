package com.lifetechnologies.solid.wt.test;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

import com.lifetechnologies.solid.wt.ReadMask;

/**
 * Unit Test for the class com.lifetechnologies.solid.wt.ReadMask
 * @author mullermw
 *
 */
public class ReadMaskTest {

	@Test
	public void test() {
		assertTrue(new ReadMask("").asBitString(50).equals            ("11111111111111111111111111111111111111111111111111"));
		assertTrue(new ReadMask("46"  ).asBitString(50).equals        ("11111111111111111111111111111111111111111111101111"));
		assertTrue(new ReadMask("46..").asBitString(50).equals        ("11111111111111111111111111111111111111111111100000"));
		assertTrue(new ReadMask("46..50").asBitString().equals        ("11111111111111111111111111111111111111111111100000"));
		assertTrue(new ReadMask("1-10").asBitString(50).equals        ("00000000001111111111111111111111111111111111111111"));            
		assertTrue(new ReadMask("1-10,15,46..").asBitString(50).equals("00000000001111011111111111111111111111111111100000"));
		assertTrue(new ReadMask("10-1,15,46..").asBitString(50).equals("00000000001111011111111111111111111111111111100000"));
		assertTrue(new ReadMask("51").asBitString(50).equals            ("11111111111111111111111111111111111111111111111111"));
	}
	
}
