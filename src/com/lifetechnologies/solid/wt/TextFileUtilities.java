package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.TreeMap;

/**
 * User: tuchbb
 * Date: Aug 22, 2008
 * Time: 8:59:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class TextFileUtilities {

    public static HashSet<String> loadSetFromFile(File fileText, String delimiter, int indexOfColumn) throws IOException {

        HashSet<String> set = new HashSet<String>();
        BufferedReader reader = new BufferedReader(new FileReader(fileText));
        String line;
        while ((line = reader.readLine()) != null) {
            String tokens[] = line.split(delimiter);
            set.add(tokens[indexOfColumn]);
        }
        reader.close();
        return set;
    }
    
    public static TreeMap<Integer, Double> loadIntegerToDoubleMapFromFile(File fileContainingData,
                                                                          int columnIndexKey, int columnIndexValue,
                                                                          String delimimiter,
                                                                          int numberOfHeaderLines) throws IOException {

        TreeMap<Integer, Double>  map = new TreeMap<Integer, Double> ();
        BufferedReader reader = new BufferedReader(new FileReader(fileContainingData));
        String line;
        for (int i = 0; i < numberOfHeaderLines; i++) reader.readLine();
        while ((line = reader.readLine()) != null) {
            String columns[] = line.split(delimimiter);
            if (columns[columnIndexValue].length() > 0) //drop missing data
                map.put(Integer.parseInt(columns[columnIndexKey]), Double.valueOf(columns[columnIndexValue]));
        }


        reader.close();
        return map;
    }
    

    public static TreeMap<String, Double> loadStringToDoubleMapFromFile(File fileContainingData,
                                                                          int columnIndexKey, int columnIndexValue,
                                                                          String delimimiter,
                                                                          int numberOfHeaderLines) throws IOException {

        TreeMap<String, Double> map = new TreeMap<String, Double>();
        BufferedReader reader = new BufferedReader(new FileReader(fileContainingData));
        String line;
        for (int i = 0; i < numberOfHeaderLines; i++) reader.readLine();
        while ((line = reader.readLine()) != null) {
            String columns[] = line.split(delimimiter);
            if (columns[columnIndexValue].length() > 0) //drop missing data
                map.put(columns[columnIndexKey], Double.valueOf(columns[columnIndexValue]));
        }


        reader.close();
        return map;
    }
    
    public static TreeMap<String, Long> loadStringToLongMapFromFile(
			File fileContainingData, int columnIndexKey, int columnIndexValue,
			String delimimiter, int numberOfHeaderLines) throws IOException {

		TreeMap<String, Long> map = new TreeMap<String, Long>();
		BufferedReader reader = new BufferedReader(new FileReader(
				fileContainingData));
		String line;
		for (int i = 0; i < numberOfHeaderLines; i++)
			reader.readLine();
		while ((line = reader.readLine()) != null) {
			String columns[] = line.split(delimimiter);
			if (columns[columnIndexValue].length() > 0) // drop missing data
				map.put(columns[columnIndexKey], Long
						.valueOf(columns[columnIndexValue]));
		}

		reader.close();
		return map;
	}

    public static TreeMap<String, String> loadStringToStringMapFromFile(File fileContainingData,
                                                                          int columnIndexKey, int columnIndexValue,
                                                                          String delimimiter,
                                                                          int numberOfHeaderLines) throws IOException {

        TreeMap<String, String> map = new TreeMap<String, String>();
        BufferedReader reader = new BufferedReader(new FileReader(fileContainingData));
        String line;
        for (int i = 0; i < numberOfHeaderLines; i++) reader.readLine();
        while ((line = reader.readLine()) != null) {
            String columns[] = line.split(delimimiter);
            if (columns[columnIndexValue].length() > 0) //drop missing data
                map.put(columns[columnIndexKey], columns[columnIndexValue]);
        }


        reader.close();
        return map;
    }

    public static String loadFirstLineOfFile(File file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        reader.close();
        return line;
    }

    public static int countNumberOfLinesInFile(File file) throws IOException {
        int count = 0;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        while (reader.readLine() != null)   count++;
        reader.close();
        return count;
    }

    public static int countNumberOfLinesInFileNotStartingWithString(File file, String prefixOfLinesToSkip) throws IOException {
        int count = 0;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null)
            if (!line.startsWith(prefixOfLinesToSkip))   count++;
        reader.close();
        return count;
    }

}
