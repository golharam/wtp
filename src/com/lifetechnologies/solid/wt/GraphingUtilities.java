package com.lifetechnologies.solid.wt;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.File;

import com.lowagie.text.DocumentException;
import com.lowagie.text.Document;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;


/**
 * User: tuchbb
 * Date: Aug 25, 2008
 * Time: 5:51:59 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class GraphingUtilities {

    /**
    * Converts a JFreeChart to PDF syntax.
    * @param filePDF	the PDF output file
    * @param chart		the JFreeChart
    * @param width		the width of the resulting PDF
    * @param height	the height of the resulting PDF
    */
    public static void convertToPdf(JFreeChart chart, int width, int height, File filePDF) {

        Document document = new Document(new com.lowagie.text.Rectangle(width, height));
        try {
            PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream(filePDF));
            document.open();
            PdfContentByte cb = writer.getDirectContent();
            PdfTemplate tp = cb.createTemplate(width, height);
            Graphics2D g2d = tp.createGraphics(width, height, new DefaultFontMapper());
            Rectangle2D r2d = new Rectangle2D.Double(0, 0, width, height);
            chart.draw(g2d, r2d);
            g2d.dispose();
            cb.addTemplate(tp, 0, 0);
        } catch(DocumentException de) {
            de.printStackTrace();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        document.close();
    }

    public static JFreeChart generateChartXYSeriesOverFusion(String idOfFusion,
                                                      String sequenceOfFusion, 
                                                      int positionOfFusionInSequence,
                                                      String nameOfFivePrimeFusionGene,
                                                      String nameOfThreePrimeFusionGene,
                                                      XYSeriesCollection dataSetForChart,
                                                      boolean zoomedImage) {

        JFreeChart chart = ChartFactory.createXYLineChart(null, idOfFusion, "Read coverage", dataSetForChart, PlotOrientation.VERTICAL, true, false, false);
        chart.getXYPlot().setBackgroundPaint(Color.WHITE);
        int lengthOfSequence = sequenceOfFusion.length();
        chart.getXYPlot().addAnnotation(new XYLineAnnotation(-20, 0, lengthOfSequence +20, 0, new BasicStroke(3), Color.BLACK));

        // add gene boxes
        double coordinatesGeneA[] = { 0,2 , positionOfFusionInSequence,2 , positionOfFusionInSequence,-2 , 0,-2 };
        double coordinatesGeneB[] = { positionOfFusionInSequence,2 , lengthOfSequence,2 , lengthOfSequence,-2, positionOfFusionInSequence,-2 };
        if (zoomedImage) {
            coordinatesGeneA = new double[] { Math.max(1,positionOfFusionInSequence-50),2 , positionOfFusionInSequence,2 , positionOfFusionInSequence,-2 , Math.max(1,positionOfFusionInSequence-50),-2 };
            coordinatesGeneB = new double[] { positionOfFusionInSequence,2 , Math.min(lengthOfSequence, positionOfFusionInSequence+50),2 , Math.min(lengthOfSequence, positionOfFusionInSequence+50),-2, positionOfFusionInSequence,-2 };
        }
        chart.getXYPlot().addAnnotation(new XYPolygonAnnotation(coordinatesGeneA, new BasicStroke(3), Color.BLACK, Color.WHITE));
        chart.getXYPlot().addAnnotation(new XYPolygonAnnotation(coordinatesGeneB, new BasicStroke(3), Color.BLACK, Color.WHITE));
        // add gene names
        XYTextAnnotation textAnnotationGeneNameA = new XYTextAnnotation(nameOfFivePrimeFusionGene, positionOfFusionInSequence /2, 0);
        XYTextAnnotation textAnnotationGeneNameB = new XYTextAnnotation(nameOfThreePrimeFusionGene, positionOfFusionInSequence + (lengthOfSequence - positionOfFusionInSequence) /2, 0);
        if (zoomedImage) {
            textAnnotationGeneNameA = new XYTextAnnotation(nameOfFivePrimeFusionGene, positionOfFusionInSequence -25, 0);
            textAnnotationGeneNameB = new XYTextAnnotation(nameOfThreePrimeFusionGene, positionOfFusionInSequence + 25, 0);
        }
        textAnnotationGeneNameA.setFont(textAnnotationGeneNameA.getFont().deriveFont(Font.BOLD, 16.0f));
        textAnnotationGeneNameA.setPaint(Color.BLUE);
        textAnnotationGeneNameB.setFont(textAnnotationGeneNameB.getFont().deriveFont(Font.BOLD, 16.0f));
        textAnnotationGeneNameB.setPaint(Color.ORANGE);
        chart.getXYPlot().addAnnotation(textAnnotationGeneNameA);
        chart.getXYPlot().addAnnotation(textAnnotationGeneNameB);
        // add gene sequences
        XYTextAnnotation textAnnotationGeneSequenceA = new XYTextAnnotation(sequenceOfFusion.substring(Math.max(0, positionOfFusionInSequence -50), positionOfFusionInSequence), positionOfFusionInSequence, -5);
        textAnnotationGeneSequenceA.setTextAnchor(TextAnchor.CENTER_RIGHT);
        textAnnotationGeneSequenceA.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 28));
        textAnnotationGeneSequenceA.setPaint(Color.BLUE);
        chart.getXYPlot().addAnnotation(textAnnotationGeneSequenceA);
        XYTextAnnotation textAnnotationGeneSequenceB = new XYTextAnnotation(sequenceOfFusion.substring(positionOfFusionInSequence, Math.min(sequenceOfFusion.length(), positionOfFusionInSequence +50)), positionOfFusionInSequence, -5);
        textAnnotationGeneSequenceB.setTextAnchor(TextAnchor.CENTER_LEFT);
        textAnnotationGeneSequenceB.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 28));
        textAnnotationGeneSequenceB.setPaint(Color.ORANGE);
        chart.getXYPlot().addAnnotation(textAnnotationGeneSequenceB);


        chart.getXYPlot().getRangeAxis().setRange(-10, 100);
        if (zoomedImage)    chart.getXYPlot().getDomainAxis().setRange(positionOfFusionInSequence - 50, positionOfFusionInSequence +50);
        else    chart.getXYPlot().getDomainAxis().setRange(-20, lengthOfSequence +20);

        chart.getXYPlot().getRangeAxis().setLabelFont(chart.getXYPlot().getRangeAxis().getLabelFont().deriveFont(Font.BOLD, 16.0f));
        chart.getXYPlot().getDomainAxis().setLabelFont(chart.getXYPlot().getDomainAxis().getLabelFont().deriveFont(Font.BOLD, 16.0f));

        chart.getXYPlot().getRangeAxis().setTickLabelFont(chart.getXYPlot().getRangeAxis().getTickLabelFont().deriveFont(Font.BOLD, 16.0f));
        chart.getXYPlot().getDomainAxis().setTickLabelFont(chart.getXYPlot().getDomainAxis().getTickLabelFont().deriveFont(Font.BOLD, 16.0f));

        TickUnits tickUnits = new TickUnits();
        tickUnits.add(new NumberTickUnit(10));
        chart.getXYPlot().getRangeAxis().setStandardTickUnits(tickUnits);

        chart.getXYPlot().setRangeGridlinesVisible(false);
        chart.getXYPlot().setDomainGridlinesVisible(false);

        chart.getXYPlot().getRenderer().setBaseStroke(new BasicStroke(3));

        //chart.getXYPlot().getRenderer().setSeriesPaint(1, Color.GREEN);

        chart.getLegend().setItemFont(chart.getLegend().getItemFont().deriveFont(Font.BOLD, 16.0f));

        chart.setBackgroundPaint(Color.WHITE);
        return chart;
    }

    public static JFreeChart generateChartXYSeriesOverSequence(String idOfSequence,
                                                               int lengthOfSequence,
                                                               XYSeriesCollection dataSetForChart) {

        JFreeChart chart = ChartFactory.createXYLineChart(null, idOfSequence, "Read coverage", dataSetForChart, PlotOrientation.VERTICAL, true, false, false);
        chart.getXYPlot().setBackgroundPaint(Color.WHITE);
        chart.getXYPlot().addAnnotation(new XYLineAnnotation(-20, 0, lengthOfSequence +20, 0, new BasicStroke(3), Color.BLACK));

        // add sequence annotation box
        double coordinatesGene[] = { 0,2 , lengthOfSequence-1,2 , lengthOfSequence-1,-2 , 0,-2 };
        chart.getXYPlot().addAnnotation(new XYPolygonAnnotation(coordinatesGene, new BasicStroke(3), Color.BLACK, Color.WHITE));
        // add sequence id
        XYTextAnnotation textAnnotationSequenceId = new XYTextAnnotation(idOfSequence, lengthOfSequence /2, 0);
        textAnnotationSequenceId.setFont(textAnnotationSequenceId.getFont().deriveFont(Font.BOLD, 16.0f));
        textAnnotationSequenceId.setPaint(Color.BLUE);
        chart.getXYPlot().addAnnotation(textAnnotationSequenceId);

        chart.getXYPlot().getRangeAxis().setRange(-10, 100);
        chart.getXYPlot().getDomainAxis().setRange(-20, lengthOfSequence +20);

        chart.getXYPlot().getRangeAxis().setLabelFont(chart.getXYPlot().getRangeAxis().getLabelFont().deriveFont(Font.BOLD, 16.0f));
        chart.getXYPlot().getDomainAxis().setLabelFont(chart.getXYPlot().getDomainAxis().getLabelFont().deriveFont(Font.BOLD, 16.0f));

        chart.getXYPlot().getRangeAxis().setTickLabelFont(chart.getXYPlot().getRangeAxis().getTickLabelFont().deriveFont(Font.BOLD, 16.0f));
        chart.getXYPlot().getDomainAxis().setTickLabelFont(chart.getXYPlot().getDomainAxis().getTickLabelFont().deriveFont(Font.BOLD, 16.0f));

        TickUnits tickUnits = new TickUnits();
        tickUnits.add(new NumberTickUnit(10));
        chart.getXYPlot().getRangeAxis().setStandardTickUnits(tickUnits);

        chart.getXYPlot().setRangeGridlinesVisible(false);
        chart.getXYPlot().setDomainGridlinesVisible(false);

        chart.getXYPlot().getRenderer().setBaseStroke(new BasicStroke(3));

        //chart.getXYPlot().getRenderer().setSeriesPaint(1, Color.GREEN);

        chart.getLegend().setItemFont(chart.getLegend().getItemFont().deriveFont(Font.BOLD, 16.0f));

        chart.setBackgroundPaint(Color.WHITE);
        return chart;
    }

}
