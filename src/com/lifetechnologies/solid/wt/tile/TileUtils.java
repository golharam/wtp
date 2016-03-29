package com.lifetechnologies.solid.wt.tile;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import com.lifetechnologies.solid.wt.Strand;

public class TileUtils {

	private TileUtils() {}
	
	public static <T extends Tile> Set<Set<T>> clusterTiles(Collection<T> tiles) {
		List<T> tileList = new ArrayList<T>(tiles);
		Set<Set<T>> clusters = new HashSet<Set<T>>();
		while (!tileList.isEmpty()) {
			T seed = tileList.remove(0);
			Set<T> cluster = new HashSet<T>();
			cluster.add(seed);
			//Recursively search for overlapping tiles.
			Stack<T> recursionStack = new Stack<T>();
			recursionStack.add(seed);
			while (!recursionStack.isEmpty()) {
				Tile tile = recursionStack.pop();
				for (Iterator<T> it = tileList.iterator(); it.hasNext(); ) {
					T candidate = it.next();
					if (tile.overlapsWith(candidate)) {
						cluster.add(candidate);
						recursionStack.push(candidate);
						it.remove();
					}
				}
			}
			clusters.add(cluster);
		}
		return clusters;
	}
	
	public static void main(String[] args) {
		Set<Tile> tiles = new HashSet<Tile>();
		tiles.add(SimpleTile.newTile("contig1", Strand.POSITIVE, 100, 150));
		tiles.add(SimpleTile.newTile("contig1", Strand.POSITIVE, 145, 200));
		tiles.add(SimpleTile.newTile("contig1", Strand.POSITIVE, 201, 300));
		tiles.add(SimpleTile.newTile("contig1", Strand.POSITIVE, 98, 101));
		tiles.add(SimpleTile.newTile("contig1", Strand.NEGATIVE, 100, 150));
		tiles.add(SimpleTile.newTile("contig2", Strand.POSITIVE, 100, 150));
		System.out.println(clusterTiles(tiles));
	}
}
