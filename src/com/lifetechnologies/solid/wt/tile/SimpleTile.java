package com.lifetechnologies.solid.wt.tile;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.lifetechnologies.solid.wt.Strand;
import com.lifetechnologies.util.EmbeddedIntegerComparator;

public class SimpleTile extends AbstractTile {

	private String contigId;
	private Strand strand;
	private int start;
	private int end;
	
	protected SimpleTile(String contigId, Strand strand, int start , int end) {
		if (contigId == null) throw new IllegalArgumentException("contigId cannot be null.");
		if (strand == null) throw new IllegalArgumentException("strand cannot be null.");
		if (start < 0) throw new IllegalArgumentException("start cannot be less than zero.");
		if (end < 0) throw new IllegalArgumentException("end cannot be less than zero.");
		if (start > end ) {
			int tmp = start;
			start = end;
			end = tmp;
		}
		this.contigId = contigId;
		this.strand = strand;
		this.start = start;
		this.end = end;
	}

	public static Tile newTile(String contigId, Strand strand, int start , int end ) {
		return new SimpleTile(contigId, strand, start, end);
	}
	
	@Override
	public String getContigId() {
		return contigId;
	}
	
	@Override
	public Strand getStrand() {
		return strand;
	}
	
	@Override
	public int getStart() {
		return start;
	}
	
	@Override
	public int getEnd() {
		return end;
	}
	
	public static void main(String[] args) {
		Tile tile1 = SimpleTile.newTile("1", Strand.POSITIVE, 10, 20);
		Tile tile2 = SimpleTile.newTile("1", Strand.POSITIVE, 15, 30);
		System.out.println(tile1.getOverlap(tile2));
		
		List<String> strings = Arrays.asList(new String[]{"chr10", "chr20", "chr2", "chr1"});
		Collections.sort(strings, new EmbeddedIntegerComparator());
		System.out.println(strings);
	}
	
}


