package com.lifetechnologies.solid.wt;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * User: tuchbb
 * Date: Nov 19, 2008
 * Time: 9:25:22 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ContiguousGenomicFeature implements Cloneable {

    private String idOfReferenceSequence;
    private Strand strand;
    private int coordinateOfStart;
    private int coordinateOfEnd;

    private double overlapFractionFromMyPerspective;
    private double overlapFractionFromItsPerspective;

    private ContiguousGenomicFeature parent;

    private Double valueA;
    private Double valueB;
    private String type;
    private String label;
    private String additionalInfo;
    private String annotation;
    
    private static final Pattern genomicLocationPattern = Pattern.compile("(\\S+?)(?::(\\d+)(?:-(\\d+))?)?");

    public ContiguousGenomicFeature(String idOfReferenceSequence, Strand strand, int coordinateOfStart, int coordinateOfEnd) {
        this.idOfReferenceSequence = idOfReferenceSequence;
        this.strand = strand;
        this.coordinateOfStart = coordinateOfStart;
        this.coordinateOfEnd = coordinateOfEnd;
    } 

    public ContiguousGenomicFeature(int coordinateOfStart, int coordinateOfEnd) {
        this.coordinateOfStart = coordinateOfStart;
        this.coordinateOfEnd = coordinateOfEnd;
    }

    protected Object clone() throws CloneNotSupportedException {
        return super.clone();    
    }

    /**
     * Output in GFF format.
     * @return
     */
    public String toString() {

        String stringIdOfReferenceSequence = ".";
        if (this.idOfReferenceSequence != null)
            stringIdOfReferenceSequence = this.idOfReferenceSequence;

        String stringLabel = ".";
        if (this.label != null)
            stringLabel = this.label;

        String stringType = ".";
        if (this.type != null)
            stringType = this.type;

        String stringValueA = ".";
        if (this.valueA != null)
            stringValueA = "" + this.valueA;

        String stringValueB = ".";
        if (this.valueB != null)
            stringValueB = "" + this.valueB;

        String stringAdditionalInfo = ".";
        if (this.additionalInfo != null)
            stringAdditionalInfo = this.additionalInfo;

        return stringIdOfReferenceSequence + "\t" + stringLabel + "\t" + stringType + "\t"
                + this.coordinateOfStart + "\t" + this.coordinateOfEnd
                + "\t" + stringValueA + "\t" + this.strand.toChar() + "\t" + stringValueB + "\t" + stringAdditionalInfo;
    }

    public String getIdOfReferenceSequence() {
        return idOfReferenceSequence;
    }

    public void setIdOfReferenceSequence(String idOfReferenceSequence) {
        this.idOfReferenceSequence = idOfReferenceSequence;
    }

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    public int getCoordinateOfStart() {
        return coordinateOfStart;
    }

    public void setCoordinateOfStart(int coordinateOfStart) {
        this.coordinateOfStart = coordinateOfStart;
    }

    public int getCoordinateOfEnd() {
        return coordinateOfEnd;
    }

    public void setCoordinateOfEnd(int coordinateOfEnd) {
        this.coordinateOfEnd = coordinateOfEnd;
    }

    public double getOverlapFractionFromMyPerspective() {
        return overlapFractionFromMyPerspective;
    }

    public void setOverlapFractionFromMyPerspective(double overlapFractionFromMyPerspective) {
        this.overlapFractionFromMyPerspective = overlapFractionFromMyPerspective;
    }

    public double getOverlapFractionFromItsPerspective() {
        return overlapFractionFromItsPerspective;
    }

    public void setOverlapFractionFromItsPerspective(double overlapFractionFromItsPerspective) {
        this.overlapFractionFromItsPerspective = overlapFractionFromItsPerspective;
    }

    public int getLength() {
        return this.coordinateOfEnd - this.coordinateOfStart +1;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getAdditionalInfo() {
        return additionalInfo;
    }

    public void setAdditionalInfo(String additionalInfo) {
        this.additionalInfo = additionalInfo;
    }

    public Double getValueA() {
        return valueA;
    }

    public void setValueA(Double valueA) {
        this.valueA = valueA;
    }

    public String getLabel() {
        return label;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public Double getValueB() {
        return valueB;
    }

    public void setValueB(Double valueB) {
        this.valueB = valueB;
    }

    public ContiguousGenomicFeature getParent() {
        return parent;
    }

    public void setParent(ContiguousGenomicFeature parent) {
        this.parent = parent;
    }

    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    public String getAnnotation() {
        return annotation;
    }

	/**
	 * Parse a string representation of features into a Set of ContiguousGenomicFeature objects.
	 * String Examples:  chr1:1-10000
	 *                   chr5:100-10000, chr10:10000-4000000
	 * @param featureString the string to parse
	 * @return parsed set of features.
	 */
    public static Set<ContiguousGenomicFeature> parseFeatures(String featureString) {
    	Set<ContiguousGenomicFeature> features = new HashSet<ContiguousGenomicFeature>();
    	if (Utilities.isBlank(featureString)) return features;
		for (String s : featureString.split("\\s*,\\s*")) {
			Matcher matcher = genomicLocationPattern.matcher(s);
			if (matcher.matches()) {
				String header = matcher.group(1);
				ContiguousGenomicFeature feature = new ContiguousGenomicFeature(header,Strand.EITHER, -1, -1);        			

				int start = -1;
				if (matcher.group(2) != null)
					start = Integer.parseInt(matcher.group(2));
				feature.setCoordinateOfStart(start);

				int end = -1;
				if (matcher.group(3) != null)
					end = Integer.parseInt(matcher.group(3));
				feature.setCoordinateOfEnd(end);

				if (feature.getCoordinateOfStart() > feature.getCoordinateOfEnd()) {
					int tmp = feature.getCoordinateOfStart();
					feature.setCoordinateOfStart(feature.getCoordinateOfEnd());
					feature.setCoordinateOfEnd(tmp);
				}
				features.add(feature);
			}
		}
		return features;
    }
    
    /**
     * 
     * @return a string representation of this feature.  This string can be
     *         parsed with parseFeatures().
     */
    public String toBriefString() {
    	StringBuffer buf = new StringBuffer(this.getIdOfReferenceSequence());
    	if (this.getCoordinateOfStart() > -1) {
    		buf.append(":");
    		buf.append(this.getCoordinateOfStart());
    		buf.append("-");
    	}
    	if (this.getCoordinateOfEnd() > -1)
    		buf.append(this.getCoordinateOfEnd());
    	return buf.toString();
    }

}
