package com.lifetechnologies.solid.wt.tile;

import com.lifetechnologies.solid.wt.Strand;

public interface Tile extends Comparable<Tile> {

	public String getContigId();
	public Strand getStrand();
	public int getStart();
	public int getEnd();
	public int getLength();
	public boolean overlapsWith(Tile t);
	public int getOverlap(Tile t);
	public int distanceFrom(Tile t);
	public boolean contains(int i);
	
}
