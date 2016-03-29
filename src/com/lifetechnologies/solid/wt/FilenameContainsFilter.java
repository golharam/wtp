package com.lifetechnologies.solid.wt;

import java.io.FilenameFilter;
import java.io.File;

/**
 * User: tuchbb
 * Date: Aug 20, 2008
 * Time: 1:37:16 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class FilenameContainsFilter implements FilenameFilter {

    String stringThatFileNameMustContain = "";

    public FilenameContainsFilter(String stringThatFileNameMustContain) {
        this.stringThatFileNameMustContain = stringThatFileNameMustContain;
    }

    public boolean accept(File dir, String name) {
        if (name.indexOf(stringThatFileNameMustContain) > -1)
            return true;
        else
            return false;
    }
}
