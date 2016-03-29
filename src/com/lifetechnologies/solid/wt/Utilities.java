package com.lifetechnologies.solid.wt;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.HashMap;

import com.lifetechnologies.util.ProcessId;

/**
 * User: tuchbb
 * Date: Aug 22, 2008
 * Time: 8:59:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class Utilities {

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
    
    public static HashMap<String, String> loadSequenceDBFromFastaFile(File fastaFile) throws IOException {

        HashMap<String, String> mapIdToSequence = new HashMap<String, String>();
        BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
        String line;
        String currentHeader = "";
        String currentSequence = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (currentHeader.length() > 0)
                    mapIdToSequence.put(currentHeader, currentSequence);
                currentHeader = line.substring(1);
                currentSequence = "";
            } else
                currentSequence += line;
        }
        if (currentHeader.length() > 0)
            mapIdToSequence.put(currentHeader, currentSequence);

        reader.close();
        return mapIdToSequence;
    }

    public static void concatenateFiles(File[] files, File fileConcatenated, boolean deleteSourceFilesWhenDone) throws IOException {

        byte[] buffer = new byte[1024 * 64];
        FileOutputStream fileOutputStreamConcatenated = new FileOutputStream(fileConcatenated);
        for (int i = 0; i < files.length; i++) {
            FileInputStream fileInputStream = new FileInputStream(files[i]);
            int numberOfBytesRead = 0;
            while ((numberOfBytesRead = fileInputStream.read(buffer)) != -1)
                fileOutputStreamConcatenated.write(buffer, 0, numberOfBytesRead);
            fileInputStream.close();
            if (deleteSourceFilesWhenDone)
                files[i].delete();
        }

        fileOutputStreamConcatenated.close();
    }
    
    /**
     * true if s is 'true' or 'yes' or 't' ignoring case.
     * @param s
     * @return
     */
    public static Boolean isTrue(String s) {
    	if (s == null) return false;
    	if (Boolean.parseBoolean(s.trim()))   return true;
    	if (s.trim().equalsIgnoreCase("yes")) return true;
    	if (s.trim().equalsIgnoreCase("t"))   return true;
    	return false;
    }

    /**
     * true if s is 'false' or 'no' or 'f', ignoring case.
     * @param s
     * @return
     */
    public static Boolean isFalse(String s) {
    	if (s == null) return false;
    	if (!Boolean.parseBoolean(s.trim())) return true;
    	if (s.trim().equalsIgnoreCase("no")) return true;
    	if (s.trim().equalsIgnoreCase("f"))  return true;
    	return false;
    }
    
    /**
     * Concatenate the toString() values of objs with delimiter.
     * @param <T>
     * @param objs
     * @param delimiter
     * @return
     */
    public static <T> String join(final Collection<T> objs, final String delimiter) {
      if (objs == null || objs.isEmpty())
        return "";
      Iterator<T> iter = objs.iterator();
      StringBuffer buffer = new StringBuffer(iter.next().toString());
      while (iter.hasNext()) {
    	T nextValue = iter.next();
    	if (nextValue == null) {
    		buffer.append(delimiter).append("null");
    	} else {
    		buffer.append(delimiter).append(nextValue.toString());
    	}
      }
      return buffer.toString();
    }
         
    /**
     * Concatenate the toString() values of objects with delimiter.
     * @param delimiter
     * @param objects
     * @return
     */
    public static String join (final String delimiter, final Object ... objects) {
    	if (objects.length == 0) return "";
    	StringBuffer buffer = new StringBuffer(objects[0].toString());
    	for (int i=1; i<objects.length; i++) {
    		Object value = objects[i];
    		if (value == null) value = "null";
    		buffer.append(delimiter).append(value);
    	}
    	return buffer.toString();
    }
    
    public static String join (final String delimiter, final String ... strings) {
    	return join(delimiter, (Object[])strings);
    }
    
    /**
     * true if s is null or all whitespace.
     * @param s
     * @return
     */
    public static boolean isBlank(final String s) {
    	if (s == null ) return true;
    	return s.trim().equals("");
    }
    
    /**
     * Print information about memory usage to the specified stream.
     * @param out
     */
    public static void printMemoryUsage(PrintStream out) {
    	out.printf("free=%12d total=%12d used=%12d \n", Runtime.getRuntime().freeMemory(), Runtime.getRuntime().totalMemory(), Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
    }

    /**
     * Convert each String in 'in' to an Integer.  Will throw a RuntimeException
     * if any of the Strings cannot be converted.
     * @param in list of Strings to convert.
     * @return list of converted Integers.
     */
    public static List<Integer> toIntegerList(Collection<String> in) {
    	List<Integer> list = new ArrayList<Integer>();
    	for (String str : in)
    		list.add(Integer.valueOf(str));
    	return list;
    }
    
    /**
     * Convert each String in 'in' to an Integer.  Will throw a RuntimeException
     * if any of the Strings cannot be converted.
     * @param in list of Strings to convert.
     * @return list of converted Integers.
     */
    public static List<Double> toDoubleList(Collection<String> in) {
    	List<Double> list = new ArrayList<Double>();
    	for (String str : in)
    		list.add(Double.valueOf(str));
    	return list;
    }
    
    /**
     * Copy an input stream to an output stream.
     * @param in
     * @param out
     * @param closeIn close the input stream when finished.
     * @param closeOut close the output stream when finished.
     * @throws IOException if an error occurs while copying.
     */
    public static void pump(InputStream in, OutputStream out, boolean closeIn, boolean closeOut) throws IOException{
        byte[] bytes = new byte[1024];
        int read;
        try{
            while((read=in.read(bytes))!= -1)
                out.write(bytes, 0, read);
        }finally{
            if(closeIn)
                in.close();
            if(closeOut)
                out.close();
        }
    }
    
    
    /**
     * Extract the stack trace of an Exception as a String.
     * @param e the Exception to convert.
     * @return the Stack trace as a string.
     */
    public static String toStackTrace(Exception e) {
    	if (e == null) return "";
    	ByteArrayOutputStream bstream = new ByteArrayOutputStream();
		PrintWriter writer = new PrintWriter(bstream);
		e.printStackTrace(writer);
		writer.close();
		return bstream.toString();
    }
    
    /**
     * Recursively delete a non-empty directory
     * 
     * from http://www.rgagnon.com/javadetails/java-0483.html
     * 
     * @param path to delete
     * @return true if directory was deleted.
     */
    static public boolean deleteDirectory(File path) {
    	if( path.exists() ) {
    		File[] files = path.listFiles();
    		for(int i=0; i<files.length; i++) {
    			if(files[i].isDirectory()) {
    				deleteDirectory(files[i]);
    			}
    			else {
    				files[i].delete();
    			}
    		}
    	}
    	return( path.delete() );
    }
    
    /**
     * Return a file that is guaranteed to be unique in this JVM
     * and avoids collisions with names returned by the same method in other 
     * JVMs by never returning a file that exists and incorporating the 
     * host name, process id and system time in the filename.
     * 
     * @param parent directory in which to create the file.
     * @param prefix prefix of the filename.
     * @return a unique file name in parent.
     */
    public static File getUniqueFileName(File parent, String prefix) {
    	if (prefix==null) prefix = "";
    	String hostName = "unknown_host";
    	try {
    		hostName = java.net.InetAddress.getLocalHost().getHostName();
    	} catch (UnknownHostException e) {}
    	File file = null;
    	while (file == null ) {
    		long now = System.currentTimeMillis();
    		String fileName = prefix.concat(hostName+"."+ProcessId.getProcessId())+"."+now;
    		file = new File(parent, fileName);
    		if (file.exists()) file = null;
   			try { Thread.sleep(1); } catch (InterruptedException e) {}
    	}
    	return file;
    }
    
    /**
     * 
     * @param <T>
     * @param collection
     * @return null if collection is null or empty.  Otherwise the first value
     * in the collection.
     */
    public static <T> T firstValue(Collection<T> collection) {
    	if (collection == null) return null;
    	Iterator<T> it = collection.iterator();
    	if (it == null) return null;
    	return it.next();
    }
    
    public static String toStringSafeForFile(String string) {
    	return string.replaceAll("[^\\w\\d._~#-]","_");
    }
    
}
