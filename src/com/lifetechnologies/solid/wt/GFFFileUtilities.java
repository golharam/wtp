package com.lifetechnologies.solid.wt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.lifetechnologies.solid.wt.splice.AnnotatedFeature;
import com.lifetechnologies.util.EmbeddedIntegerComparator;



/**
 * User: tuchbb
 * Date: Dec 1, 2008
 * Time: 2:16:06 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class GFFFileUtilities {

	//private static Logger logger = Logger.getLogger(GFFFileUtilities.class.toString());
	
    /**
     * Note that coordinates in GFF files are expected to use a 1-based system.
     * 
     * @param fileFeatureAnnotationGFF
     * @param nameOfChromosome
     * @param strandToLoadFrom
     * @param numberOfBasesToLengthenFeaturesBy
     * @return
     * @throws IOException
     */
    public static LinkedHashMap<Integer, ContiguousGenomicFeature[]> loadCoordinateToFeaturesMapFromGFF(File fileFeatureAnnotationGFF,
                                                                                                        String nameOfChromosome,
                                                                                                        Strand strandToLoadFrom,
                                                                                                        int numberOfBasesToLengthenFeaturesBy) throws Exception {

        LinkedHashMap<Integer, ContiguousGenomicFeature[]> mapCoordinateToFeature = new LinkedHashMap<Integer, ContiguousGenomicFeature[]>();

        BufferedReader readerFaetureAnnotationFile = new BufferedReader(new FileReader(fileFeatureAnnotationGFF));
        String line;
        while ((line = readerFaetureAnnotationFile.readLine()) != null) {
            if (!line.startsWith("#") && !line.startsWith("track ")) {
                String tokens[] = line.split("\t");

                Strand strand = Strand.toStrand(tokens[6].toCharArray()[0]);

                if (tokens[0].equalsIgnoreCase(nameOfChromosome) && strand == strandToLoadFrom) {

                    int coordinateStart = Integer.parseInt(tokens[3]);
                    coordinateStart = Math.max(1, coordinateStart - numberOfBasesToLengthenFeaturesBy);
                    int coordinateEnd = Integer.parseInt(tokens[4]);
                    coordinateEnd += numberOfBasesToLengthenFeaturesBy;

                    //ContiguousGenomicFeature feature = new ContiguousGenomicFeature(coordinateStart, coordinateEnd);
                    ContiguousGenomicFeature feature = new ContiguousGenomicFeature(nameOfChromosome, strand, coordinateStart, coordinateEnd);       // too big?
                    feature.setLabel(tokens[1]);
                    feature.setType(tokens[2]);
                    feature.setAdditionalInfo(tokens[8]);
                    for (int currentCoordinate = coordinateStart; currentCoordinate <= coordinateEnd; currentCoordinate++) {
                        if (mapCoordinateToFeature.containsKey(currentCoordinate)) {
                            ContiguousGenomicFeature features[] = mapCoordinateToFeature.get(currentCoordinate);
                            ContiguousGenomicFeature featuresNew[] = Arrays.copyOf(features, features.length +1);
                            featuresNew[features.length] = feature;
                            mapCoordinateToFeature.put(currentCoordinate, featuresNew);
                        } else {
                            ContiguousGenomicFeature featuresNew[] = new ContiguousGenomicFeature[1];
                            featuresNew[0] = feature;
                            mapCoordinateToFeature.put(currentCoordinate, featuresNew);
                        }

                    }
                }
            }
        }
        readerFaetureAnnotationFile.close();

        return mapCoordinateToFeature;

    }
    
    /**
     * Parse a GFF stream.
     * Note that coordinates in GFF files are expected to use a 1-based system.
     * 
     * @param fileFeatureAnnotationGFF the stream to parse
     * @return list of features.
     * @throws IOException
     */
    public static List<ContiguousGenomicFeature> parseGff(InputStream fileFeatureAnnotationGFF) throws IOException, ParseException {

        List<ContiguousGenomicFeature> featureList = new ArrayList<ContiguousGenomicFeature>();
        
        BufferedReader readerFaetureAnnotationFile = new BufferedReader(new InputStreamReader(fileFeatureAnnotationGFF));
        try {
	        String line;
	        int lineCount = 0;
	        while ((line = readerFaetureAnnotationFile.readLine()) != null) {
	        	lineCount++;
	            if (!line.startsWith("#") && !line.startsWith("track ")) {
	                String tokens[] = line.split("\t");
	
	                Strand strand = null;
	                try {
	                	strand = Strand.toStrand(tokens[6].toCharArray()[0]);
	                } catch (Exception e) {
	                	throw new ParseException("Unable to parse strand on line# " + lineCount + ":\n" + line, lineCount );
	                }
	                String nameOfChromosome = tokens[0];
	
	                int coordinateStart = Integer.parseInt(tokens[3]);
	                int coordinateEnd = Integer.parseInt(tokens[4]);
	
	                ContiguousGenomicFeature feature = new ContiguousGenomicFeature(nameOfChromosome, strand, coordinateStart, coordinateEnd);       // too big?
	                feature.setLabel(tokens[1]);
	                feature.setType(tokens[2]);
	                feature.setAdditionalInfo(tokens[8]);
	                
	                featureList.add(feature);
	                
	            }
	        }
        } finally {
        	readerFaetureAnnotationFile.close();
        }
        return featureList;

    }
    
    /**
     * Returns a sorted list of exons parsed from inputStream
     * */
    public static List<AnnotatedFeature> parseUniqueExonsFromGtf(InputStream inputStream) throws IOException, ParseException {
    	if (inputStream == null) throw new IllegalArgumentException("inputStream cannot be null.");
        Map<AnnotatedFeature, AnnotatedFeature> exons = new HashMap<AnnotatedFeature, AnnotatedFeature>();
        BufferedReader reader = null;
        try {
        	reader = new BufferedReader(new InputStreamReader(inputStream));
        	int lineCount = 0;
        	for (String line = reader.readLine(); line != null; line = reader.readLine()) {
        		lineCount++;
	            if (line.startsWith("#") || line.startsWith("track ")) continue;
	            String tokens[] = line.split("\t");
	            if (tokens.length < 8) {
	            	throw new ParseException("Unable to parse line# " + lineCount + ":\n" + line, lineCount );
	            }
                Strand strand = null;
                try {
                	strand = Strand.toStrand(tokens[6].toCharArray()[0]);
                } catch (Exception e) {
                	throw new ParseException("Unable to parse strand on line# " + lineCount + ":\n" + line, lineCount );
                }
                String nameOfChromosome = tokens[0];
                String type = tokens[2];
                if (!type.equals("exon")) continue;
                int coordinateStart = Integer.parseInt(tokens[3]);
                int coordinateEnd = Integer.parseInt(tokens[4]);

                AnnotatedFeature prototype = new AnnotatedFeature(nameOfChromosome, strand, coordinateStart, coordinateEnd);
                	
                prototype.setLabel(tokens[1]);
                prototype.setType(tokens[2]);

                AnnotatedFeature feature = exons.get(prototype);
                if (feature == null) {
                	feature = prototype;
                	exons.put(prototype, feature);
                }
                
                GtfAttributes attributes = GtfAttributes.parseGtfAttributes(tokens[8]);
                feature.addTranscript(attributes.geneId, attributes.transcriptId);
        	}
        } finally {
        	if (reader != null) reader.close();
        }
        List<AnnotatedFeature> exonList= new ArrayList<AnnotatedFeature>();
        exonList.addAll(exons.keySet());
        Collections.sort(exonList, new Comparator<AnnotatedFeature>() {
	    	@Override
	    	public int compare(AnnotatedFeature o1, AnnotatedFeature o2) {
	    		if (!o1.getIdOfReferenceSequence().equals(o2.getIdOfReferenceSequence())) 
	    			return EmbeddedIntegerComparator.INSTANCE.compare(o1.getIdOfReferenceSequence(), o2.getIdOfReferenceSequence());
	    		if (o1.getCoordinateOfStart() != o2.getCoordinateOfStart())
	    			return o1.getCoordinateOfStart() - o2.getCoordinateOfStart();
	    		if (o1.getCoordinateOfEnd() != o2.getCoordinateOfEnd())
	    			return o1.getCoordinateOfEnd() - o2.getCoordinateOfEnd();
	    		if (o1.getStrand() != o2.getStrand())
	    			return o1.getStrand().compareTo(o2.getStrand());
	    		return 0;
	    	}
        });
        return exonList;
    }

	public static void main(String[] args) throws Exception {
		File gtfFile = new File("//siena/mullermw/test/1.2_testing/refGene.my.gtf");
		FileInputStream stream = new FileInputStream(gtfFile);
		for (AnnotatedFeature feature : parseUniqueExonsFromGtf(stream)) {
			System.out.println(feature);
		}
	}
	
    public static void transformExonGFFIntoGeneGFF(File fileExonGFF, File fileGeneGFF) throws Exception {

        HashMap<String, ContiguousGenomicFeature> mapGeneNameToGene = new HashMap<String, ContiguousGenomicFeature>();
        HashSet<String> setOfGeneNamesOnMultipleContigs = new HashSet<String>();
        HashSet<String> setOfGeneNamesOnBothStrands = new HashSet<String>();
        BufferedReader readerExonGFF = new BufferedReader(new FileReader(fileExonGFF));
        String line;
        while ((line = readerExonGFF.readLine()) != null) {
            if (line.length() > 0 && !line.startsWith("#")) {
                String tokens[] = line.split("\t");
                String nameOfContig = tokens[0];
                String nameOfGene = tokens[8].substring(3, tokens[8].indexOf(';'));
                int coordinateStart = Integer.parseInt(tokens[3]);
                int coordinateEnd = Integer.parseInt(tokens[4]);

                Strand strand = Strand.toStrand(tokens[6].toCharArray()[0]);

                if (mapGeneNameToGene.containsKey(nameOfGene)) {
                    ContiguousGenomicFeature gene = mapGeneNameToGene.get(nameOfGene);
                    if (!gene.getIdOfReferenceSequence().equals(nameOfContig)) {
                        setOfGeneNamesOnMultipleContigs.add(nameOfGene);
                        System.out.println("ERROR: " + nameOfGene + " apparently has exons on two different chromosomes.");
                    } else if (gene.getStrand() != strand) {
                        setOfGeneNamesOnBothStrands.add(nameOfGene);
                        System.out.println("ERROR: " + nameOfGene + " apparently has exons on both strands.");
                    } else {
                        if (gene.getCoordinateOfStart() > coordinateStart)
                            gene.setCoordinateOfStart(coordinateStart);
                        if (gene.getCoordinateOfEnd() < coordinateEnd)
                            gene.setCoordinateOfEnd(coordinateEnd);
                    }
                } else {
                    ContiguousGenomicFeature gene = new ContiguousGenomicFeature(nameOfContig, strand, coordinateStart, coordinateEnd);
                    gene.setAdditionalInfo("ID=" + nameOfGene);
                    gene.setLabel(tokens[1] + "2Gene");
                    gene.setType("gene");
                    mapGeneNameToGene.put(nameOfGene, gene);
                }
            }
        }
        readerExonGFF.close();

        System.out.println(setOfGeneNamesOnMultipleContigs.size() + " genes have exons on multiple contigs.");
        System.out.println(setOfGeneNamesOnBothStrands.size() + " genes have exons on both strands.");

        BufferedWriter writerGeneGFF = new BufferedWriter(new FileWriter(fileGeneGFF));
        Iterator<String> iteratorOverGeneNames = mapGeneNameToGene.keySet().iterator();
        while (iteratorOverGeneNames.hasNext()) {
            String nameOfGene = iteratorOverGeneNames.next();
            writerGeneGFF.write(mapGeneNameToGene.get(nameOfGene).toString());
            writerGeneGFF.newLine();
        }
        writerGeneGFF.close();
    }
    

	
	static class GtfAttributes {
		public final String geneId;
		public final String transcriptId;
		public static final Pattern geneIdAttributePattern = Pattern.compile("gene_id\\s*\"\\s*([^\"\\s]+?)\\s*\"\\s*;");
		public static final Pattern transcriptIdAttributePattern = Pattern.compile(geneIdAttributePattern.pattern().replaceFirst("gene_id", "transcript_id"));
		
		
		public GtfAttributes(String geneId, String transcriptId) {
			if (geneId == null) throw new IllegalArgumentException("geneId cannot be null.");
			if (transcriptId == null) throw new IllegalArgumentException("transcriptId cannot be null");
			this.geneId = geneId;
			this.transcriptId = transcriptId;
		}

		public static GtfAttributes parseGtfAttributes(String attributesField) {
			if (attributesField == null) throw new IllegalArgumentException("attributesField cannot be null.");
			String geneId, transcriptId;
			Matcher matcher = geneIdAttributePattern.matcher(attributesField);
			if (!matcher.find()) throw new IllegalArgumentException("Cannot parse '"+attributesField+"' as GTF attributes.");
			geneId = matcher.group(1);
			matcher = transcriptIdAttributePattern.matcher(attributesField);
			if (!matcher.find()) throw new IllegalArgumentException("Cannot parse '"+attributesField+"' as GTF attributes.");
			transcriptId = matcher.group(1);
			return new GtfAttributes(geneId, transcriptId);
		}
		
	}
}

