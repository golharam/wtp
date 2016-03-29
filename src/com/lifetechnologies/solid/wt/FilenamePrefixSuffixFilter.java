package com.lifetechnologies.solid.wt;

import java.io.File;
import java.io.FilenameFilter;

/**
 * Created by IntelliJ IDEA.
 * Revision: $Rev$
 * User: Tuch
 */
public class FilenamePrefixSuffixFilter implements FilenameFilter {

    private String prefix;
    private String suffix;

    public FilenamePrefixSuffixFilter(String prefix, String suffix) {
        this.prefix = prefix;
        this.suffix = suffix;
    }

    public boolean accept(File dir, String name) {
        return (name.startsWith(prefix) && name.endsWith(suffix));
    }
}
