package com.lifetechnologies.solid.wt;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.plot.PlotOrientation;

import java.util.*;
import java.io.File;
import java.io.Serializable;
import java.text.DecimalFormat;

/**
 * User: tuchbb
 * Date: Oct 12, 2008
 * Time: 9:35:49 PM
 * Revision: $Rev$
 *
 */
public class Histogram implements Serializable {

	public static final long serialVersionUID = 1;
	
    private int numberOfBins;
    private double sizeOfBins;
    private double valueOfFirstBin;
    private String nameOfBinSeries;

    private double bins[];
    private HashMap<String, double[]> mapHistogramNameToFrequencies = new HashMap<String, double[]>();

    boolean reportRelativeFrequencies = false;

    public Histogram(int numberOfBins, double sizeOfBins, double valueOfFirstBin, String nameOfBinSeries) {
        this.numberOfBins = numberOfBins;
        this.sizeOfBins = sizeOfBins;
        this.valueOfFirstBin = valueOfFirstBin;
        this.nameOfBinSeries = nameOfBinSeries;

        this.bins = new double[numberOfBins];
        for (int i = 0; i < numberOfBins; i++)
            this.bins[i] = valueOfFirstBin + i * sizeOfBins;

    }

    public void addHistogramSeries(String name, double[] histogramSeries) {
        this.mapHistogramNameToFrequencies.put(name, histogramSeries);
    }

    public double[] getHistogramSeries(String name) {
        return this.mapHistogramNameToFrequencies.get(name);
    }

    public Set<String> getHistogramSeriesNames() {
        return this.mapHistogramNameToFrequencies.keySet();
    }

    public int getNumberOfBins() {
        return this.numberOfBins;
    }

    public double getSizeOfBins() {
        return this.sizeOfBins;
    }

    public double getValueOfFirstBin() {
        return this.valueOfFirstBin;
    }

    public double[] getBins() {
        return this.bins;
    }

    public String getNameOfBinSeries() {
        return nameOfBinSeries;
    }

    public void setNameOfBinSeries(String nameOfBinSeries) {
        this.nameOfBinSeries = nameOfBinSeries;
    }

    public String toString() {

        String histogramsAsString = "";

        String namesOfHistogramSeries[] = new String[this.getHistogramSeriesNames().size()];
        namesOfHistogramSeries = this.getHistogramSeriesNames().toArray(namesOfHistogramSeries);
        Arrays.sort(namesOfHistogramSeries);

        double totalCountsForEachSeries[] = new double[namesOfHistogramSeries.length];
        for (int i = 0; i < totalCountsForEachSeries.length; i++) {
            totalCountsForEachSeries[i] = 0;
            for (int j = 0; j < this.getHistogramSeries(namesOfHistogramSeries[i]).length; j++)
                totalCountsForEachSeries[i] += this.getHistogramSeries(namesOfHistogramSeries[i])[j];
        }

        if (this.reportRelativeFrequencies)
            histogramsAsString += "RELATIVE FREQUENCIES\n";
        else
            histogramsAsString += "ABSOLUTE FREQUENCIES\n";

        histogramsAsString += this.nameOfBinSeries.replaceAll(" ", "_");
        for (int i = 0; i < namesOfHistogramSeries.length; i++)
            histogramsAsString += "\t" + (namesOfHistogramSeries[i] + " (N = " + (long)totalCountsForEachSeries[i] + ")").replaceAll(" ", "_");
        histogramsAsString += "\n";

        for (int i = 0; i < this.bins.length; i++) {
            histogramsAsString += this.bins[i];
            for (int j = 0; j < namesOfHistogramSeries.length; j++) {       // (new DecimalFormat("##0.0")).format(100.0 * countsOfBestReadMappingByScore[i] / countTotalBestReadMappings) + "%\t"
                if (reportRelativeFrequencies)  histogramsAsString +=  "\t" + (new DecimalFormat("##0.0")).format(100.0 * (this.getHistogramSeries(namesOfHistogramSeries[j])[i] / totalCountsForEachSeries[j])) + "%";
                else    histogramsAsString +=  "\t" + (long)this.getHistogramSeries(namesOfHistogramSeries[j])[i];
            }
            histogramsAsString += "\n";
        }

        histogramsAsString += "Total";
        if (this.reportRelativeFrequencies) {
            for (int i = 0; i < totalCountsForEachSeries.length; i++)
                histogramsAsString += "\t100.0%";
        } else {
            for (int i = 0; i < totalCountsForEachSeries.length; i++)
                histogramsAsString += "\t" + (long)totalCountsForEachSeries[i];
        }
        histogramsAsString += "\n";

        return histogramsAsString;
    }

    public void toPDF(File fileOutputImage, boolean centerBins, boolean cumulative, boolean logScaleYAxis, Double fractionOfYAxisToDisplayForTailFocus) throws Exception {

        String namesOfHistogramSeries[] = new String[this.getHistogramSeriesNames().size()];
        namesOfHistogramSeries = this.getHistogramSeriesNames().toArray(namesOfHistogramSeries);
        Arrays.sort(namesOfHistogramSeries);
        XYSeriesCollection data = new XYSeriesCollection();

        double totalCountsForEachSeries[] = new double[namesOfHistogramSeries.length];
        for (int i = 0; i < totalCountsForEachSeries.length; i++) {
            totalCountsForEachSeries[i] = 0;
            for (int j = 0; j < this.getHistogramSeries(namesOfHistogramSeries[i]).length; j++)
                totalCountsForEachSeries[i] += this.getHistogramSeries(namesOfHistogramSeries[i])[j];
            if (logScaleYAxis)
                totalCountsForEachSeries[i] += this.getHistogramSeries(namesOfHistogramSeries[i]).length;
        }

        for (int i = 0; i < namesOfHistogramSeries.length; i++) {
            double sum = 0;
            XYSeries series = new XYSeries(namesOfHistogramSeries[i] + " (N = " + (long)totalCountsForEachSeries[i] + ")");
            for (int j = 0; j < bins.length; j++) {
                double scoreBin = j * this.sizeOfBins + this.valueOfFirstBin;                
                if (centerBins) // center each value in range of values spanned by the bin
                    scoreBin -= this.sizeOfBins / 2.0;  // FIXME: this was previosuly adding half the bin size, but should be subtracting this value, right?
                    //scoreBin += this.sizeOfBins / 2.0;

                double frequencyOfCurrentBin = this.getHistogramSeries(namesOfHistogramSeries[i])[j];
                if (logScaleYAxis)
                    frequencyOfCurrentBin++;
                if (this.reportRelativeFrequencies)
                    frequencyOfCurrentBin /= totalCountsForEachSeries[i];
                sum += frequencyOfCurrentBin;

                if (cumulative)
                    series.add(scoreBin, sum);
                else 
                    series.add(scoreBin, frequencyOfCurrentBin);
            }
            data.addSeries(series);
        }

        String labelYAxis = "Frequency";
        if (this.reportRelativeFrequencies)  labelYAxis = "Relative Frequency";

        JFreeChart chart = ChartFactory.createXYLineChart(null, this.nameOfBinSeries, labelYAxis, data, PlotOrientation.VERTICAL, true, false, false);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, true);
        chart.getXYPlot().setRenderer(renderer);

        chart.setBackgroundPaint(chart.getXYPlot().getBackgroundPaint());

        if (logScaleYAxis)
            chart.getXYPlot().setRangeAxis(new LogarithmicAxis(labelYAxis));

        if (fractionOfYAxisToDisplayForTailFocus != null && fractionOfYAxisToDisplayForTailFocus < 1.0 && this.reportRelativeFrequencies) {    // FIXME: add focusOnTail functionaility when this.reportRelativeFrequencies = false
            if (cumulative)
                chart.getXYPlot().getRangeAxis().setLowerBound(1.0 - fractionOfYAxisToDisplayForTailFocus);
            else
                chart.getXYPlot().getRangeAxis().setUpperBound(fractionOfYAxisToDisplayForTailFocus);
        }

        //ChartUtilities.saveChartAsPNG(fileOutputImage, chart, 600, 400);
        GraphingUtilities.convertToPdf(chart, 600, 400, fileOutputImage);
    }

    public void createNewHistogramSeries(String nameOfHistogramSeries) {
        this.mapHistogramNameToFrequencies.put(nameOfHistogramSeries, new double[this.bins.length]);
    }

    public void addObservationToHistogramSeries(String nameOfHistogramSeries, double valueObserved) {
        int indexOfBin = (int)Math.ceil((valueObserved - this.valueOfFirstBin) / this.sizeOfBins);
        indexOfBin = Math.max(indexOfBin, 0);
        indexOfBin = Math.min(indexOfBin, bins.length -1);
        if (!this.mapHistogramNameToFrequencies.containsKey(nameOfHistogramSeries))
            this.mapHistogramNameToFrequencies.put(nameOfHistogramSeries, new double[this.numberOfBins]);
        this.mapHistogramNameToFrequencies.get(nameOfHistogramSeries)[indexOfBin]++;
    }

    public boolean isReportRelativeFrequencies() {
        return this.reportRelativeFrequencies;
    }

    public void setReportRelativeFrequencies(boolean reportRelativeFrequencies) {
        this.reportRelativeFrequencies = reportRelativeFrequencies;
    }

    public void addAllSeriesFromHistogram(Histogram histogram) {

        Set<String> namesOfSeries = histogram.getHistogramSeriesNames();
        Iterator<String> iteratorNamesOfSeries = namesOfSeries.iterator();
        while (iteratorNamesOfSeries.hasNext()) {
            String nameOfSeries = iteratorNamesOfSeries.next();
            if (this.getHistogramSeries(nameOfSeries) == null)
                this.addHistogramSeries(nameOfSeries, histogram.getHistogramSeries(nameOfSeries));
            else {
                double series[] = this.getHistogramSeries(nameOfSeries);
                double seriesToAdd[] = histogram.getHistogramSeries(nameOfSeries);
                for (int i = 0; i < series.length; i++)
                    series[i] += seriesToAdd[i];
                                
            }
        }
    }

    public static void main(String args[]) throws Exception {
        Histogram histogramOfAllelicRatios = new Histogram(21, 1.0, -10.0, "Allelic Ratios [ log2(ref/SNP) ]");
        String nameOfSeriesAlleleRatiosSampleA = "Allelic Ratios Sample A";
        histogramOfAllelicRatios.createNewHistogramSeries(nameOfSeriesAlleleRatiosSampleA);
        histogramOfAllelicRatios.addObservationToHistogramSeries(nameOfSeriesAlleleRatiosSampleA, 0.1);
        histogramOfAllelicRatios.addObservationToHistogramSeries(nameOfSeriesAlleleRatiosSampleA, -1.2);
        histogramOfAllelicRatios.addObservationToHistogramSeries(nameOfSeriesAlleleRatiosSampleA, -12);
        histogramOfAllelicRatios.addObservationToHistogramSeries(nameOfSeriesAlleleRatiosSampleA, 15.0);
        histogramOfAllelicRatios.setReportRelativeFrequencies(true);
        histogramOfAllelicRatios.toPDF(new File("/home/tuchbb/temp.pdf"), true, false, false, null);
    }
}

