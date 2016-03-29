package com.lifetechnologies.solid.wt.mapper;

import java.io.PrintWriter;
import java.text.FieldPosition;
import java.text.Format;
import java.text.ParsePosition;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.Strand;
import com.lifetechnologies.solid.wt.Utilities;

public class GffFeature extends ContiguousGenomicFeature implements Comparable<GffFeature> {
	private int score;
	private int phase = -1;
	private String featureType;
	private String source;
	private String beadId;
	private int indexOfReference;
	
	public static final Format FORMAT = new Format() {
		public static final long serialVersionUID = 1;
		
		@Override
		public StringBuffer format(Object obj, StringBuffer toAppendTo,
				FieldPosition pos) {
			if (obj instanceof GffFeature == false) throw new IllegalArgumentException("Not a GffFeature");
			return toAppendTo.append(((GffFeature)obj).asGffRow());
		}
		
		@Override
		public Object parseObject(String source, ParsePosition pos) {
			return new GffFeature(source.substring(pos.getIndex()));
		}
		
		@Override 
		public Object parseObject(String source) {
			return new GffFeature(source);
		}
	};
	
	public GffFeature(String idOfReference, Strand strand, int coordinateOfStart, int coordinateOfEnd) {
		super(idOfReference == null  ? "" : idOfReference, strand == null ? Strand.EITHER : strand, coordinateOfStart, coordinateOfEnd);
		if (coordinateOfStart < 1) throw new IllegalArgumentException("coordinateOfStart cannot be less than 1");
		if (coordinateOfEnd < coordinateOfStart) throw new IllegalArgumentException("coordinateOfEnd cannot be less than coordinateOfStart");
	}
	
	public GffFeature(String gffLine) {
		super(null, null, 1, 1);
		String[] fields = gffLine.split("\\s*\t\\s*");
		this.setIdOfReferenceSequence(fields[0]);
		this.setSource(fields[1]);
		this.setFeatureType(fields[2]);
		this.setCoordinateOfStart(Integer.parseInt(fields[3]));
		this.setCoordinateOfEnd(Integer.parseInt(fields[4]));
		this.setScore(Integer.parseInt(fields[5]));
		this.setStrand(Strand.toStrand(fields[6].charAt(0)));
		if (!fields[7].equals(".")) {
			this.setPhase(Integer.parseInt(fields[7]));
		}
		this.setAdditionalInfo(fields[8]);
	}

	public int getScore() {
		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public int getPhase() {
		return phase;
	}

	public void setPhase(int phase) {
		this.phase = phase;
	}

	public String getFeatureType() {
		return featureType;
	}

	public void setFeatureType(String featureType) {
		this.featureType = featureType;
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}
	
	@Override
	public void setAdditionalInfo(String additionalInfo) {
		super.setAdditionalInfo(additionalInfo);
		String[] attributes = additionalInfo.split(";");
		this.beadId = attributes[0].substring(3);
		this.indexOfReference = Integer.parseInt(attributes[4].substring(2));
	}
	
	public void toGffRow(PrintWriter out) {
		if (out == null) return;
		out.println(this.asGffRow());
	}
	
	public String asGffRow() {
		return Utilities.join("\t", 
			this.getIdOfReferenceSequence() == null ? "null" : this.getIdOfReferenceSequence(), 
			this.getSource(), 
			this.getFeatureType(), 
			this.getCoordinateOfStart(), 
			this.getCoordinateOfEnd(), 
			this.getScore(), 
			this.getStrand().toChar(), 
			this.getPhase() > -1 ? this.getPhase() : ".", 
			this.getAdditionalInfo());
	}
	
	@Override
	public int compareTo(GffFeature that) {
		if (this == that ) return 0;
		if (this.equals(that)) return 0;
		if (this.indexOfReference != that.indexOfReference) return this.indexOfReference - that.indexOfReference;
		if (this.getCoordinateOfStart() != that.getCoordinateOfStart()) return this.getCoordinateOfStart() - that.getCoordinateOfStart();
		if (this.getCoordinateOfEnd() != that.getCoordinateOfEnd()) return this.getCoordinateOfEnd() - that.getCoordinateOfEnd();
		if (this.getStrand() != that.getStrand()) return this.getStrand().compareTo(that.getStrand());
		return this.beadId.compareTo(that.beadId);
	}
	
}

