package com.lifetechnologies.solid.wt;

public class GenomicPosition implements Comparable<GenomicPosition> {

	private int sequenceIndex;
	private int sequencePosition;
	private Strand strand;
	
	public GenomicPosition(int sequenceIndex, int sequencePosition, Strand strand) {
		if (sequenceIndex < 1) throw new IllegalArgumentException("sequenceIndex cannot be less than 1.");
		if (sequencePosition < 0) throw new IllegalArgumentException("sequencePosition cannot be less than 0.");
		if (strand == null) strand = Strand.EITHER;
		this.sequenceIndex = sequenceIndex;
		this.sequencePosition = sequencePosition;
		this.strand = strand;
	}
	
	public GenomicPosition(int sequenceIndex, int sequencePosition) {
		this(sequenceIndex, 
			 sequencePosition < 0 ? Math.abs(sequencePosition) : sequencePosition, 
			 sequencePosition < 0 ? Strand.NEGATIVE : Strand.POSITIVE);
	}
	
	public GenomicPosition derivePosition(int distance) {
		int newPos = this.sequencePosition + distance;
		if (newPos < 0 ) throw new IllegalArgumentException("Illegal Attempt to derive a position less than zero.");
		return new GenomicPosition(this.sequenceIndex, newPos, this.strand);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) return false;
		if (obj == this) return true;
		if (! (obj instanceof GenomicPosition)) return false;
		GenomicPosition that = (GenomicPosition)obj;
		return this.sequenceIndex == that.sequenceIndex &&
			   this.sequencePosition == that.sequencePosition &&
			   this.strand == that.strand;
	}
	
	@Override
	public String toString() {
		return GenomicPosition.class.toString() + ":"+this.sequenceIndex+","+this.sequencePosition+","+this.strand;
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
	
	@Override
	public int compareTo(GenomicPosition that) {
		if (that == null ) return -1;
		if (this.sequenceIndex != that.sequenceIndex) return this.sequenceIndex - that.sequenceIndex;
		if (this.sequencePosition != that.sequencePosition) return this.sequencePosition - that.sequencePosition;
		return this.strand.compareTo(that.strand);
	}

	public int getSequenceIndex() {
		return sequenceIndex;
	}

	public int getSequencePosition() {
		return sequencePosition;
	}

	public Strand getStrand() {
		return strand;
	}
}
