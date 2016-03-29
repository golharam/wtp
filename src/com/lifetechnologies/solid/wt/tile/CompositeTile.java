package com.lifetechnologies.solid.wt.tile;

import java.util.SortedSet;
import java.util.TreeSet;

import com.lifetechnologies.solid.wt.Strand;

public class CompositeTile extends AbstractTile {

	private SortedSet<Tile> tiles = new TreeSet<Tile>();

	public void addTile(Tile tile) {
		if (tile == null) throw new IllegalArgumentException("tile cannot be null.");
		if (tiles.contains(tile)) return;
		if (tiles.size() > 0 && !tile.getContigId().equals(tiles.first().getContigId()))
			throw new IllegalArgumentException("Cannot add Tiles from different contigs.");
		if (tiles.size() > 0 && !tile.getStrand().equals(tiles.first().getStrand()))
			throw new IllegalArgumentException("Cannot add Tiles from different strands.");
		tiles.add(tile);
	}
	
	@Override
	public String getContigId() {
		if (tiles.size() == 0) return null;
		return tiles.first().getContigId();
	}
	
	@Override
	public Strand getStrand() {
		if (tiles.size() == 0) return null;
		return tiles.first().getStrand();
	}
	
	@Override
	public int getStart() {
		return tiles.first().getStart();
	}
	
	@Override
	public int getEnd() {
		return tiles.last().getEnd();
	}
	
}
