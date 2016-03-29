package com.lifetechnologies.solid.wt.test;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Map;

import org.junit.Test;

import com.lifetechnologies.util.FileHeader;

public class FileHeaderTest {
	
	@Test
	public void testParseHeader() throws IOException {
		StringBuffer buf = new StringBuffer();
		buf.append("#foo=bar\n");
		buf.append("data\n");
		buf.append("#bar=foo\n");
		buf.append("#empty=\n");
		buf.append("#int=42\n");
		buf.append("#double=0.42\n");
		Map<String, FileHeader.Value> map = FileHeader.parseHeader(new ByteArrayInputStream(buf.toString().getBytes()));
		assertEquals(new FileHeader.Value("bar"), map.get("foo"));
		assertEquals(new FileHeader.Value("foo"), map.get("bar"));
		assertEquals(new FileHeader.Value(""), map.get("empty"));
		assertFalse(map.get("foo").isNumeric());
		assertTrue(map.get("int").isNumeric());
		assertEquals((long)42, (long)map.get("int").intValue());
		assertEquals(0.42, map.get("double").doubleValue(),0.001);
		FileHeader.writeHeader(map, System.out, "");
	}
}
