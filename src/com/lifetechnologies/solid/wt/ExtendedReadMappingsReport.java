package com.lifetechnologies.solid.wt;

import java.text.DecimalFormat;
import java.io.*;

import com.lifetechnologies.solid.wt.mapper.FilteringMode;

/**
 * User: tuchbb
 * Date: Oct 10, 2008
 * Time: 6:16:10 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExtendedReadMappingsReport implements Serializable {

	public static final long serialVersionUID = 1;
	
    private Integer minMappingsRequiredBeforeReportingRead;
    private Integer maxMappingsAllowedBeforeNotReportingRead;
    private Integer minAlignmentScoreForReportingMapping;
    private Integer minScoreGapToSecondBestAlignmentForUniqueness;

    private Integer countOfReadsProcessed;
    private Integer countOfReadsWithAtLeastOneSplitHavingAFilterTag;
    private Integer countOfReadsWithAllSplitsHavingAFilterTag;

    private Integer countOfReadsWithTooFewMappings;
    private Integer countOfReadsWithTooManyMappings;
    private Integer countOfReadsWithNumberOfMappingsInRequiredRange;

    private Integer countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore;

    private Integer countOfReadsUniquelyAlignedWithMinAlignScore;

    private int[] countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit;

    private Integer countOfReadsWithReadSplitsAligningToTheSameLocation;

    private Histogram histogramOfMappingsByScore;
    
    private FilteringMode filteringMode = FilteringMode.ONE_OR_MORE;

	public String toString() {
    	ByteArrayOutputStream baos = new ByteArrayOutputStream();
    	PrintWriter writer = new PrintWriter(baos);
    	writeReport(writer);
    	writer.close();
    	return baos.toString();
    }
    
    public void writeReport(PrintWriter out) {
    	out.println("--------------\nALIGNMENT REPORT\n--------------");
    	out.println();
    	out.println("Counts:");
    	out.printf("Reads mapped:\t%d  (100%%)\n", this.getTotalMappedReads() );
   		out.printf("Reads filtered:\t%d\t(%.1f%%)",
    		this.getCountReadsFiltered(), 
    		(double)this.getCountReadsFiltered() / getTotalMappedReads() * 100);
   		out.println();
   		out.printf("Reads with too many mappings (N > %d):\t%d\t(%.1f%%)",
   				this.maxMappingsAllowedBeforeNotReportingRead,
   				this.countOfReadsWithTooManyMappings,
   				(double)this.countOfReadsWithTooManyMappings / getTotalMappedReads() * 100);
   		out.println();
   		out.printf("Reads with number of mappings in proper range (N <= %d):\t%d\t(%.1f%%)",
   				this.maxMappingsAllowedBeforeNotReportingRead,
   				this.countOfReadsWithNumberOfMappingsInRequiredRange,
   				(double)this.countOfReadsWithNumberOfMappingsInRequiredRange / getTotalMappedReads() * 100);
   		out.println();
   		out.printf("Reads with number of mappings in proper range (N <= %d) and score >= %d :\t%d\t(%.1f%%)",
   					this.maxMappingsAllowedBeforeNotReportingRead,
   					this.minAlignmentScoreForReportingMapping,
   					this.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore,
   					(double)this.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore / getTotalMappedReads() * 100);
   		out.println();
   		out.printf("Reads uniquely aligned (minScoreGapToSecondBestAlignment = %d ) with align score >= %d :\t%d\t(%.1f%%)",
   				this.minScoreGapToSecondBestAlignmentForUniqueness,
   				this.minAlignmentScoreForReportingMapping,
   				this.countOfReadsUniquelyAlignedWithMinAlignScore,
   				(double)this.countOfReadsUniquelyAlignedWithMinAlignScore / getTotalMappedReads() * 100);
   		out.println();
   		
        if (this.histogramOfMappingsByScore != null) {
            this.histogramOfMappingsByScore.setReportRelativeFrequencies(false);
            out.println();
            out.println(this.histogramOfMappingsByScore.toString());
            out.println();
            this.histogramOfMappingsByScore.setReportRelativeFrequencies(true);
            out.println(this.histogramOfMappingsByScore.toString());
            out.println();
        }
   		
    }

    public String getReport1_1() { 
        String report = "--------------\nALIGNMENT REPORT\n--------------\n\n";

        if (this.countOfReadsProcessed != null) {
            report += "Count of reads processed:\t" + this.countOfReadsProcessed + "\n";

            if (this.minMappingsRequiredBeforeReportingRead != null && this.countOfReadsWithAtLeastOneSplitHavingAFilterTag != null)
                report += "Count of reads with at least one read split filtered:\t" + this.countOfReadsWithAtLeastOneSplitHavingAFilterTag
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithAtLeastOneSplitHavingAFilterTag / this.countOfReadsProcessed) + "%)\n";

            if (this.minMappingsRequiredBeforeReportingRead != null && this.countOfReadsWithAllSplitsHavingAFilterTag != null)
                report += "Count of reads with all splits filtered:\t" + this.countOfReadsWithAllSplitsHavingAFilterTag
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithAllSplitsHavingAFilterTag / this.countOfReadsProcessed) + "%)\n";

            if (this.minMappingsRequiredBeforeReportingRead != null && this.countOfReadsWithTooFewMappings != null)
                report += "Count of reads with too few mappings (N < " + this.minMappingsRequiredBeforeReportingRead + "):\t" + this.countOfReadsWithTooFewMappings
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithTooFewMappings / this.countOfReadsProcessed) + "%)\n";

            if (this.maxMappingsAllowedBeforeNotReportingRead != null && this.countOfReadsWithTooManyMappings != null)
                report += "Count of reads with too many mappings (N > " + this.maxMappingsAllowedBeforeNotReportingRead + "):\t" + this.countOfReadsWithTooManyMappings
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithTooManyMappings / this.countOfReadsProcessed) + "%)\n";

            if (this.minMappingsRequiredBeforeReportingRead != null && this.maxMappingsAllowedBeforeNotReportingRead != null && this.countOfReadsWithNumberOfMappingsInRequiredRange != null)
                report += "Count of reads with number of mappings in proper range (" +  this.minMappingsRequiredBeforeReportingRead +  " <= N <= " + this.maxMappingsAllowedBeforeNotReportingRead + "):\t" + this.countOfReadsWithNumberOfMappingsInRequiredRange
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithNumberOfMappingsInRequiredRange / this.countOfReadsProcessed) + "%)\n";

            //report += "\n";
            if (this.minMappingsRequiredBeforeReportingRead != null && this.maxMappingsAllowedBeforeNotReportingRead != null && this.countOfReadsWithNumberOfMappingsInRequiredRange != null)
                report += "Count of reads with number of mappings in proper range (" +  this.minMappingsRequiredBeforeReportingRead +  " <= N <= " + this.maxMappingsAllowedBeforeNotReportingRead
                            + ") and align score >= " + this.minAlignmentScoreForReportingMapping + ":\t" + this.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore / this.countOfReadsProcessed) + "%)\n";

            if (this.countOfReadsUniquelyAlignedWithMinAlignScore != null && this.minScoreGapToSecondBestAlignmentForUniqueness != null)
                report += "Count of reads uniquely aligned (minScoreGapToSecondBestAlignment = " +  this.minScoreGapToSecondBestAlignmentForUniqueness
                            + ") with align score >= " + this.minAlignmentScoreForReportingMapping + ":\t" + this.countOfReadsUniquelyAlignedWithMinAlignScore
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsUniquelyAlignedWithMinAlignScore / this.countOfReadsProcessed) + "%)\n";


            if (this.countOfReadsWithReadSplitsAligningToTheSameLocation != null)
                report += "Count of aligned reads (min scoring in proper range) with read splits aligning to the same location:\t" + this.countOfReadsWithReadSplitsAligningToTheSameLocation
                            + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countOfReadsWithReadSplitsAligningToTheSameLocation / this.countOfReadsProcessed) + "%)\n";
                        
            if (this.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit != null)
                for (int i = 0; i < countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit.length; i++)
                    report += "Count of aligned reads (min scoring in proper range) for read split " + i + ":\t" + this.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[i]
                                + "\t(" + (new DecimalFormat("##0.0")).format(100.0 * this.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[i] / this.countOfReadsProcessed) + "%)\n";


        }

        if (this.histogramOfMappingsByScore != null) {
            this.histogramOfMappingsByScore.setReportRelativeFrequencies(false);
            report += "\n" + this.histogramOfMappingsByScore.toString();
            report += "\n";
            this.histogramOfMappingsByScore.setReportRelativeFrequencies(true);
            report += this.histogramOfMappingsByScore.toString();
            report += "\n";
        }

        return report;
    }


    public Integer getCountOfReadsProcessed() {
        return countOfReadsProcessed;
    }

    public void setCountOfReadsProcessed(Integer countOfReadsProcessed) {
        this.countOfReadsProcessed = countOfReadsProcessed;
    }

    public Integer getMinMappingsRequiredBeforeReportingRead() {
        return minMappingsRequiredBeforeReportingRead;
    }

    public void setMinMappingsRequiredBeforeReportingRead(Integer minMappingsRequiredBeforeReportingRead) {
        this.minMappingsRequiredBeforeReportingRead = minMappingsRequiredBeforeReportingRead;
    }

    public Integer getMaxMappingsAllowedBeforeNotReportingRead() {
        return maxMappingsAllowedBeforeNotReportingRead;
    }

    public void setMaxMappingsAllowedBeforeNotReportingRead(Integer maxMappingsAllowedBeforeNotReportingRead) {
        this.maxMappingsAllowedBeforeNotReportingRead = maxMappingsAllowedBeforeNotReportingRead;
    }

    public Integer getCountOfReadsWithNumberOfMappingsInRequiredRange() {
        return countOfReadsWithNumberOfMappingsInRequiredRange;
    }

    public void setCountOfReadsWithNumberOfMappingsInRequiredRange(Integer countOfReadsWithNumberOfMappingsInRequiredRange) {
        this.countOfReadsWithNumberOfMappingsInRequiredRange = countOfReadsWithNumberOfMappingsInRequiredRange;
    }

    public Integer getCountOfReadsWithTooFewMappings() {
        return countOfReadsWithTooFewMappings;
    }

    public void setCountOfReadsWithTooFewMappings(Integer countOfReadsWithTooFewMappings) {
        this.countOfReadsWithTooFewMappings = countOfReadsWithTooFewMappings;
    }

    public Integer getCountOfReadsWithTooManyMappings() {
        return countOfReadsWithTooManyMappings;
    }

    public void setCountOfReadsWithTooManyMappings(Integer countOfReadsWithTooManyMappings) {
        this.countOfReadsWithTooManyMappings = countOfReadsWithTooManyMappings;
    }

    public Integer getCountOfReadsWithReadSplitsAligningToTheSameLocation() {
        return countOfReadsWithReadSplitsAligningToTheSameLocation;
    }

    public void setCountOfReadsWithReadSplitsAligningToTheSameLocation(Integer countOfReadsWithReadSplitsAligningToTheSameLocation) {
        this.countOfReadsWithReadSplitsAligningToTheSameLocation = countOfReadsWithReadSplitsAligningToTheSameLocation;
    }

    public void setCountOfReadsWithAllSplitsHavingAFilterTag(Integer countOfReadsWithAllSplitsHavingAFilterTag) {
        this.countOfReadsWithAllSplitsHavingAFilterTag = countOfReadsWithAllSplitsHavingAFilterTag;
    }

    public Integer getCountOfReadsWithAllSplitsHavingAFilterTag() {
        return countOfReadsWithAllSplitsHavingAFilterTag;
    }

    public Integer getCountOfReadsWithAtLeastOneSplitHavingAFilterTag() {
        return countOfReadsWithAtLeastOneSplitHavingAFilterTag;
    }

    public void setCountOfReadsWithAtLeastOneSplitHavingAFilterTag(Integer countOfReadsWithAtLeastOneSplitHavingAFilterTag) {
        this.countOfReadsWithAtLeastOneSplitHavingAFilterTag = countOfReadsWithAtLeastOneSplitHavingAFilterTag;
    }

    public Histogram getHistogramOfMappingsByScore() {
        return histogramOfMappingsByScore;
    }

    public void setHistogramOfMappingsByScore(Histogram histogramOfMappingsByScore) {
        this.histogramOfMappingsByScore = histogramOfMappingsByScore;
    }

    public Integer getMinAlignmentScoreForReportingMapping() {
        return minAlignmentScoreForReportingMapping;
    }

    public void setMinAlignmentScoreForReportingMapping(Integer minAlignmentScoreForReportingMapping) {
        this.minAlignmentScoreForReportingMapping = minAlignmentScoreForReportingMapping;
    }

    public Integer getCountOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore() {
        return countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore;
    }

    public void setCountOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore(Integer countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore) {
        this.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore = countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore;
    }

    public Integer getCountOfReadsUniquelyAlignedWithMinAlignScore() {
        return countOfReadsUniquelyAlignedWithMinAlignScore;
    }

    public void setCountOfReadsUniquelyAlignedWithMinAlignScore(Integer countOfReadsUniquelyAlignedWithMinAlignScore) {
        this.countOfReadsUniquelyAlignedWithMinAlignScore = countOfReadsUniquelyAlignedWithMinAlignScore;
    }

    public Integer getMinScoreGapToSecondBestAlignmentForUniqueness() {
        return minScoreGapToSecondBestAlignmentForUniqueness;
    }

    public void setMinScoreGapToSecondBestAlignmentForUniqueness(Integer minScoreGapToSecondBestAlignmentForUniqueness) {
        this.minScoreGapToSecondBestAlignmentForUniqueness = minScoreGapToSecondBestAlignmentForUniqueness;
    }
    
    public FilteringMode getFilteringMode() {
		return filteringMode;
	}

	public void setFilteringMode(FilteringMode filteringMode) {
		this.filteringMode = filteringMode;
	}
	
	public int getCountReadsFiltered() {
	
		if (filteringMode == FilteringMode.OFF)	return 0;
		if (filteringMode == FilteringMode.BOTH) return this.getCountOfReadsWithAllSplitsHavingAFilterTag();
		if (filteringMode == FilteringMode.ONE_OR_MORE) return this.getCountOfReadsWithAtLeastOneSplitHavingAFilterTag();
		return -1;
	}

    public void writeReportToTextFile(File fileMappingReportText) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(fileMappingReportText));
        String lines[] = this.toString().split("\n");
        for (int i = 0; i < lines.length; i++) {
            writer.write(lines[i]);
            writer.newLine();
        }
        writer.close();
    }

    public void setCountsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit(int[] countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit) {
        this.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit = countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit;
    }

    public int[] getCountsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit() {
        return countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit;
    }

    public static ExtendedReadMappingsReport combineSerializedMappingReports(File[] filesSerializedReportOfExtendedReadMappings) throws Exception {
        ExtendedReadMappingsReport reportsOfExtendedReadMappings[] = new ExtendedReadMappingsReport[filesSerializedReportOfExtendedReadMappings.length];
        for (int i = 0; i < filesSerializedReportOfExtendedReadMappings.length; i++) {
            ObjectInputStream objectInputStream = new ObjectInputStream(new FileInputStream(filesSerializedReportOfExtendedReadMappings[i]));
            reportsOfExtendedReadMappings[i] = (ExtendedReadMappingsReport)objectInputStream.readObject();
            objectInputStream.close();
        }


        return combineMappingReports(reportsOfExtendedReadMappings);
    }

    public static ExtendedReadMappingsReport combineMappingReports(ExtendedReadMappingsReport[] reportsOfExtendedReadMappings) throws Exception {

        ExtendedReadMappingsReport reportOfExtendedReadMappingsFull = new ExtendedReadMappingsReport();

        reportOfExtendedReadMappingsFull.minMappingsRequiredBeforeReportingRead = reportsOfExtendedReadMappings[0].minMappingsRequiredBeforeReportingRead;
        reportOfExtendedReadMappingsFull.maxMappingsAllowedBeforeNotReportingRead = reportsOfExtendedReadMappings[0].maxMappingsAllowedBeforeNotReportingRead;
        reportOfExtendedReadMappingsFull.minAlignmentScoreForReportingMapping = reportsOfExtendedReadMappings[0].minAlignmentScoreForReportingMapping;
        reportOfExtendedReadMappingsFull.minScoreGapToSecondBestAlignmentForUniqueness = reportsOfExtendedReadMappings[0].minScoreGapToSecondBestAlignmentForUniqueness;
                
        reportOfExtendedReadMappingsFull.countOfReadsProcessed = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithAtLeastOneSplitHavingAFilterTag = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithAllSplitsHavingAFilterTag = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithTooFewMappings = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithTooManyMappings = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithNumberOfMappingsInRequiredRange = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore = 0;
        reportOfExtendedReadMappingsFull.countOfReadsUniquelyAlignedWithMinAlignScore = 0;
        reportOfExtendedReadMappingsFull.countOfReadsWithReadSplitsAligningToTheSameLocation = 0;
        reportOfExtendedReadMappingsFull.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit = new int[reportsOfExtendedReadMappings[0].countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit.length];
        Histogram histogramFirst = reportsOfExtendedReadMappings[0].getHistogramOfMappingsByScore();
        reportOfExtendedReadMappingsFull.histogramOfMappingsByScore = new Histogram(histogramFirst.getNumberOfBins(), histogramFirst.getSizeOfBins(), histogramFirst.getValueOfFirstBin(), histogramFirst.getNameOfBinSeries());

        for (int i = 0; i < reportsOfExtendedReadMappings.length; i++) {
            reportOfExtendedReadMappingsFull.countOfReadsProcessed += reportsOfExtendedReadMappings[i].countOfReadsProcessed;
            reportOfExtendedReadMappingsFull.countOfReadsWithAtLeastOneSplitHavingAFilterTag += reportsOfExtendedReadMappings[i].countOfReadsWithAtLeastOneSplitHavingAFilterTag;
            reportOfExtendedReadMappingsFull.countOfReadsWithAllSplitsHavingAFilterTag += reportsOfExtendedReadMappings[i].countOfReadsWithAllSplitsHavingAFilterTag;
            reportOfExtendedReadMappingsFull.countOfReadsWithTooFewMappings += reportsOfExtendedReadMappings[i].countOfReadsWithTooFewMappings;
            reportOfExtendedReadMappingsFull.countOfReadsWithTooManyMappings += reportsOfExtendedReadMappings[i].countOfReadsWithTooManyMappings;
            reportOfExtendedReadMappingsFull.countOfReadsWithNumberOfMappingsInRequiredRange += reportsOfExtendedReadMappings[i].countOfReadsWithNumberOfMappingsInRequiredRange;
            reportOfExtendedReadMappingsFull.countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore += reportsOfExtendedReadMappings[i].countOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScore;
            reportOfExtendedReadMappingsFull.countOfReadsUniquelyAlignedWithMinAlignScore += reportsOfExtendedReadMappings[i].countOfReadsUniquelyAlignedWithMinAlignScore;
            reportOfExtendedReadMappingsFull.countOfReadsWithReadSplitsAligningToTheSameLocation += reportsOfExtendedReadMappings[i].countOfReadsWithReadSplitsAligningToTheSameLocation;

            for (int j = 0; j < reportsOfExtendedReadMappings[i].countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit.length; j++)
                reportOfExtendedReadMappingsFull.countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[j] +=
                                reportsOfExtendedReadMappings[i].countsOfReadsWithNumberOfMappingsInRequiredRangeAndMinAlignScoreBySequenceSplit[j];

            reportOfExtendedReadMappingsFull.histogramOfMappingsByScore.addAllSeriesFromHistogram(reportsOfExtendedReadMappings[i].getHistogramOfMappingsByScore());

        }
        return reportOfExtendedReadMappingsFull;
    }

    private int getTotalMappedReads() {
    	return this.countOfReadsProcessed - this.countOfReadsWithTooFewMappings;
    }
    
}


