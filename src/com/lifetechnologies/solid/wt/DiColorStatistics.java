package com.lifetechnologies.solid.wt;

import java.util.HashMap;

/**
 * User: tuchbb
 * Date: Dec 15, 2008
 * Time: 3:17:28 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class DiColorStatistics {

    HashMap<String, Long> mapDiColorToCount;
    long countOfAllDiColors;
    long countOfAllPseudoCounts;

    public DiColorStatistics(long pseudoCounts) {
        this.mapDiColorToCount = new HashMap<String, Long>();
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                this.mapDiColorToCount.put(i + "" + j, pseudoCounts);
        this.countOfAllDiColors = 0;
        this.countOfAllPseudoCounts = 16 * pseudoCounts;
    }

    public Long getCountForDiColor(String diColor) {
        if (this.mapDiColorToCount.containsKey(diColor))
            return this.mapDiColorToCount.get(diColor);
        return null;
    }

    public Long incrementCountForDiColor(String diColor, long increment) {
        this.countOfAllDiColors += increment;
        return this.mapDiColorToCount.put(diColor, this.mapDiColorToCount.get(diColor) + increment);
    }

    public long getCountOfAllDiColors() {
        return countOfAllDiColors;
    }

    public long getCountAllDiColorsPlusPseudoCounts() {
        return countOfAllDiColors + countOfAllPseudoCounts;
    }

}
