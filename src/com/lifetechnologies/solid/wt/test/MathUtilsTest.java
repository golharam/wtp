package com.lifetechnologies.solid.wt.test;

import static org.junit.Assert.assertEquals;

import java.text.ParseException;

import org.junit.Test;

import static com.lifetechnologies.util.MathUtils.*;
import static org.junit.Assert.*;

public class MathUtilsTest {

    @Test
    public void testToBytes() throws Exception {
    	assertEquals(100, toBytes("100"));
    	assertEquals(100, toBytes("100b"));
    	assertEquals(KILOBYTE, toBytes("1KB"));
    	assertEquals(MEGABYTE, toBytes("1MB"));
    	assertEquals(GIGABYTE, toBytes("1gB"));
    	ParseException e = null;
    	try {
    		toBytes("1GBC");
    	} catch (ParseException ex) {
    		e = ex;
    	}
    	assertNotNull(e);
    }
    
}
