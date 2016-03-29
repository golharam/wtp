package com.lifetechnologies.solid.wt.tile;

import com.lifetechnologies.util.EmbeddedIntegerComparator;

public abstract class AbstractTile implements Tile {

	@Override
	public int getLength() {
		return getEnd() - getStart() + 1;
	}
	
	@Override
	public int distanceFrom(Tile tile) {
		if (tile == null) throw new IllegalArgumentException("tile cannot be null.");
		if (this.overlapsWith(tile)) return 0;
		int[] distances = {this.getStart() - tile.getStart(), this.getStart() - tile.getEnd(), this.getEnd() - tile.getStart(), this.getEnd() - tile.getEnd()};
		int min = Integer.MAX_VALUE;
		for (int distance : distances ) {
			if (Math.abs(distance) < Math.abs(min)) min = distance;
		}
		return min;
	}
	
	@Override
	public int getOverlap(Tile tile) {
		if (!this.overlapsWith(tile)) return -1;
		if (this.getStart() <= tile.getStart()) {
			if (this.getEnd() <= tile.getEnd()) {
				return this.getEnd() - tile.getStart() + 1;
			} else {
				return tile.getLength();
			} 
		} else {
			if (this.getEnd() <= tile.getEnd()) {
				return this.getLength();
			} else {
				return tile.getEnd() - this.getStart() + 1;
			}
		}
	}
	
	@Override
	public boolean overlapsWith(Tile tile) {
		if (this.getContigId().equals(tile.getContigId()) == false) return false;
		if (this.getStrand().equals(tile.getStrand()) == false) return false;
		return (this.getStart() >= tile.getStart() && this.getStart() <= tile.getEnd()) ||
			   (this.getEnd() >= tile.getStart() && this.getEnd() <= tile.getEnd()) ||
			   (this.getStart() <= tile.getStart() && this.getEnd() >= tile.getEnd());
	}
	
	@Override
	public boolean contains(int i) {
		return i >= this.getStart() && i <= this.getEnd();
	}
	
	@Override
	public int compareTo(Tile tile) {
		if (tile == null) return 1;
		if (this == tile) return 0;
		if (!this.getContigId().equals(tile.getContigId()))
			return EmbeddedIntegerComparator.INSTANCE.compare(this.getContigId(), tile.getContigId());
		if (this.getStart() != tile.getStart()) return this.getStart() - tile.getStart();
		if (this.getEnd() != tile.getEnd()) return this.getEnd() - tile.getEnd();
		return this.getStrand().compareTo(tile.getStrand());
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName()+":"+getContigId()+"/"+getStart()+"-"+getEnd()+"("+getStrand().toChar()+")";
	}
	
	@Override
	public int hashCode() {
		return this.toString().hashCode();
	}

}
