package com.lifetechnologies.util;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

public class FileUtils {
	//http://www.exampledepot.com/egs/java.io/DeleteDir.html
    // Deletes all files and subdirectories under dir.
    // Returns true if all deletions were successful.
    // If a deletion fails, the method stops attempting to delete and returns false.
    public static boolean deleteDir(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteDir(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
    
        // The directory is now empty so delete it
        return dir.delete();
    }
    
    public static boolean deleteFiles(File dir, FilenameFilter filter) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteFiles(new File(dir, children[i]), filter);
                if (!success) {
                    return false;
                }
            }
        }
        if (dir.isDirectory() == false && filter.accept(dir.getParentFile(), dir.getName())) return dir.delete();
        return true;
    }
    
    public static void main(String[] args) throws IOException {
    	File dir = new File("dir");
    	File subDir = new File(dir, "dir");
    	subDir.mkdirs();
    	new File(dir, "myFile.txt").createNewFile();
    	new File(dir, "myFile.out").createNewFile();
    	new File(subDir, "myFile.txt").createNewFile();
    	new File(subDir, "myFile.out").createNewFile();
    	System.out.println(deleteFiles(dir, new FilenameFilter() {
    		@Override
    		public boolean accept(File dir, String name) {
    			// TODO Auto-generated method stub
    			return name.endsWith(".out");
    		}
    	}));
    	
    }
}
