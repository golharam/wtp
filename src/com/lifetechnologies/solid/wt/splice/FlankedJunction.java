package com.lifetechnologies.solid.wt.splice;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.lifetechnologies.solid.wt.Strand;
import com.lifetechnologies.util.SeqUtils;

public class FlankedJunction implements Comparable<FlankedJunction> {
	
	private int leftStart, leftEnd, rightStart, rightEnd;
	private StringBuffer leftSequence = new StringBuffer();
	private StringBuffer rightSequence = new StringBuffer();
	private boolean isKnown;
	private Gene gene;
	
	public FlankedJunction(Gene gene) {
		if (gene == null) throw new IllegalArgumentException("gene cannot be null");
		if (!(gene.getStrand() == Strand.POSITIVE || gene.getStrand() == Strand.NEGATIVE)) throw new IllegalArgumentException("Strand must be defined.");
		this.gene = gene;
	}
	
	public StringBuffer getLeftSequence() {
		return leftSequence;
	}
	public StringBuffer getRightSequence() {
		return rightSequence;
	}
	public int getLeftStart() {
		return leftStart;
	}
	public void setLeftStart(int leftStart) {
		this.leftStart = leftStart;
	}
	public int getLeftEnd() {
		return leftEnd;
	}
	public void setLeftEnd(int leftEnd) {
		this.leftEnd = leftEnd;
	}
	public int getLeftLength() {
		return Math.abs(this.leftEnd - this.leftStart) + 1;
	}
	public int getRightStart() {
		return rightStart;
	}
	public void setRightStart(int rightStart) {
		if (rightStart <= this.leftEnd) throw new IllegalArgumentException(String.format("rightStart(%d) <= leftEnd(%d)", rightStart, this.leftEnd));
		this.rightStart = rightStart;
	}
	public int getRightEnd() {
		return rightEnd;
	}
	public void setRightEnd(int rightEnd) {
		this.rightEnd = rightEnd;
	}
	public int getRightLength() {
		return Math.abs(this.rightEnd - this.rightStart) + 1;
	}
	public int getLength() {
		return this.getLeftLength() + this.getRightLength();
	}
	public boolean isKnown() {
		return isKnown;
	}
	public void setKnown(boolean isKnown) {
		this.isKnown = isKnown;
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public String getChromId() {
		return gene.getChromId();
	}
	
	public Strand getStrand() {
		return gene.getStrand();
	}
	
	public String getSequence() {
		if (gene.getStrand() == Strand.NEGATIVE) {
			//Revcomp junctions on negative strand.
			String leftSequence = SeqUtils.reverseComplement(getRightSequence().toString());
			String rightSequence = SeqUtils.reverseComplement(getLeftSequence().toString());
			return leftSequence + rightSequence;
		} else {
			return getLeftSequence().toString() + getRightSequence().toString();
		}
	}
	
	/**
	 * 
	 * @param positionInJunctionReference zero-based position in junction seq
	 * @return zero-based position in source genomic reference.
	 */
	public int getPositionInGenomicReference(int positionInJunctionReference) {
		int positionInReference = -1;
		if (this.getStrand() == Strand.POSITIVE) {
			if (positionInJunctionReference < this.getLeftLength()) {
				positionInReference = this.getLeftStart() - 1 + positionInJunctionReference;
			} else {
				positionInReference = this.getRightStart() - 1 + positionInJunctionReference - this.getLeftLength();
			}
		} else {
			if (positionInJunctionReference < this.getRightLength()) {
				positionInReference = this.getRightEnd() - 1 - positionInJunctionReference;
			} else {
				positionInReference = this.getLeftEnd() - 1 - positionInJunctionReference + this.getRightLength();
			}
		}
		return positionInReference;
	}
	
	/**  @return Zero based position of the last position in the acceptor in junction
	 *  sequence.
	 *
	 */
	public int getLastPositionInDonorSequence() {
		if (this.getStrand() == Strand.POSITIVE) {
			return this.getLeftLength() - 1;
		}
		return this.getRightLength() -1;
	}
	
	@Override
	public int hashCode() {
		return 	(FlankedJunction.class+gene.toString()+(leftStart)+(leftEnd)+(rightStart)+(rightEnd)+"").hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) return false;
		
		if (!(obj instanceof FlankedJunction)) return false;
		FlankedJunction that = (FlankedJunction)obj;
		return      this.gene.equals(that.gene) &&
		       this.leftEnd == that.leftEnd &&
		       this.leftStart == that.leftStart &&
		       this.rightStart == that.rightStart &&
		       this.rightEnd == that.rightEnd;
	}
	
	@Override
	public int compareTo(FlankedJunction that) {
		if (that == null) return -1;
		if (this == that) return 0;
		if (!this.gene.getChromId().equals(that.getGene().getChromId())) return this.gene.getChromId().compareTo(that.gene.getChromId());
		if (this.getLeftStart() != that.getLeftStart()) return this.getLeftStart() - that.getLeftStart();
		if (this.getLeftEnd() != that.getLeftEnd()) return this.getLeftEnd() - that.getLeftEnd();
		if (this.getRightStart() != that.getRightStart()) return this.getRightStart() - that.getRightStart();
		if (this.getRightEnd() != that.getRightEnd()) return this.getRightEnd() - that.getRightEnd();
		return this.gene.getStrand().compareTo(that.gene.getStrand());
	}
	
	private static Pattern junctionStringPattern = Pattern.compile("junction:(.*?):(\\d+)-(\\d+)\\|(\\d+)-(\\d+):(.*?):([-+]):(putative|known)");
	/* 1=contigId, 2=leftStart, 3=leftEnd, 4=rightStart, 5=rightEnd, 6=geneId, 7=strand, 8=putativeOrKnown */
	
	public String toString() {
		return "junction:"+gene.getChromId()+":"+leftStart+"-"+leftEnd+"|"+rightStart+"-"+rightEnd+":"+gene.getId()+":"+gene.getStrand().toChar()+":"+(isKnown ? "known" : "putative");
	}
	
	public String getDescription() {
		return gene.getChromId()+":"+leftStart+"-"+leftEnd+"|"+rightStart+"-"+rightEnd+":"+":"+gene.getId()+gene.getStrand().toChar()+":"+(isKnown ? "known" : "putative");
	}
	
	public static FlankedJunction parseFlankedJunction(String s) {
		Matcher matcher = junctionStringPattern.matcher(s);
		if (!matcher.matches()) throw new IllegalArgumentException(s +" cannot be parsed as a FlankedJunction.");
		Gene gene = new Gene(matcher.group(6));
		gene.setStrand(Strand.toStrand(matcher.group(7).charAt(0)));
		gene.setChromId(matcher.group(1));
		FlankedJunction junction = new FlankedJunction(gene);
		junction.setKnown(matcher.group(8).trim().toLowerCase().equals("known"));
		junction.setLeftEnd(Integer.parseInt(matcher.group(3)));
		junction.setLeftStart(Integer.parseInt(matcher.group(2)));
		junction.setRightEnd(Integer.parseInt(matcher.group(5)));
		junction.setRightStart(Integer.parseInt(matcher.group(4)));
		return junction;
	}
	
	public static void main(String[] args) {
		Gene gene = new Gene("mygene");
		gene.setStrand(Strand.POSITIVE);
		gene.setChromId("chr1");
		FlankedJunction junction = new FlankedJunction(gene);
		junction.setKnown(true);
		junction.setLeftStart(70);
		junction.setLeftEnd(100);
		junction.setRightStart(1002);
		junction.setRightEnd(1050);
		String string = junction.toString();
		System.out.println(string);
		FlankedJunction parsedJunction = FlankedJunction.parseFlankedJunction(string);
		System.out.println(parsedJunction);
		System.out.println(junction.equals(parsedJunction));
		
			
		junction = FlankedJunction.parseFlankedJunction("junction:chr17_contig6:431876-431921|432418-432463:SNIP:+:known");
		for (int i=0; i<junction.getLength(); i++) {
			System.out.printf("%d -> %d\n", i, junction.getPositionInGenomicReference(i));
		}
		System.out.println("\n\n\n");
		junction = FlankedJunction.parseFlankedJunction("junction:chr17_contig6:431876-431921|432418-432463:SNIP:-:known");
		for (int i=0; i<junction.getLength(); i++) {
			System.out.printf("%d -> %d\n", i, junction.getPositionInGenomicReference(i));
		}
		
	}
}
