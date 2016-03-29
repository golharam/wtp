package com.lifetechnologies.solid.wt.test;

import java.io.File;

import static org.junit.Assert.*;
import org.junit.Test;

import com.lifetechnologies.solid.wt.mapper.SplitReadsTask;

/**
 * Unit tests for the class SpltReadTask
 * @author mullermw
 *
 */
public class ReadSplitterTest {
		
	@Test
	public void testReadSplitter() throws Exception {
		File readsFile = new File("test/input/HBR.chr17_6.100k.mixed.csfasta");
		File destDir = new File("ignore/test_output");
		destDir.mkdirs();
		int lengthLeft = 25;
		int lengthRight = 30;
		SplitReadsTask splitter = new SplitReadsTask(readsFile, destDir, lengthLeft, lengthRight);
		System.out.println(splitter);
		splitter.doTask();
		splitter = new SplitReadsTask(readsFile, destDir, lengthLeft, lengthRight);
		
		assertTrue(splitter.isDone());
	}
}
