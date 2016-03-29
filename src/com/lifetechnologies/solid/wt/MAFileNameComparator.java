package com.lifetechnologies.solid.wt;

import java.util.Comparator;
import java.io.File;

/**
 * User: tuchbb
 * Date: Sep 10, 2008
 * Time: 3:54:10 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class MAFileNameComparator implements Comparator<File> {

    public int compare(File fileA, File fileB) {

        if (fileA.equals(fileB))
            return 0;

        int indexOfLastPeriodInFilenameA = fileA.getName().lastIndexOf('.');
        int indexOfSecondToLastPeriodInFilenameA = fileA.getName().lastIndexOf('.', indexOfLastPeriodInFilenameA -1);
        int indexOfLastPeriodInFilenameB = fileB.getName().lastIndexOf('.');
        int indexOfSecondToLastPeriodInFilenameB = fileB.getName().lastIndexOf('.', indexOfLastPeriodInFilenameB -1);

        String prefixFilenameA = fileA.getName().substring(0, indexOfSecondToLastPeriodInFilenameA);
        String prefixFilenameB = fileB.getName().substring(0, indexOfSecondToLastPeriodInFilenameB);

        if (prefixFilenameA.equalsIgnoreCase(prefixFilenameB)) {
            String readIndicesFilenameA = fileA.getName().substring(indexOfSecondToLastPeriodInFilenameA +1, indexOfLastPeriodInFilenameA);
            String readIndicesFilenameB = fileB.getName().substring(indexOfSecondToLastPeriodInFilenameB +1, indexOfLastPeriodInFilenameB);
            int indexOfFirstReadA = Integer.parseInt(readIndicesFilenameA.split("_")[0]);
            int indexOfFirstReadB = Integer.parseInt(readIndicesFilenameB.split("_")[0]);

            if (indexOfFirstReadA > indexOfFirstReadB)  return 1;
            else    return -1;

        } else
            return prefixFilenameA.compareTo(prefixFilenameB);


    }
}
