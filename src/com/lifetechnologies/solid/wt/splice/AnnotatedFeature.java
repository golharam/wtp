package com.lifetechnologies.solid.wt.splice;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.Strand;

public class AnnotatedFeature extends ContiguousGenomicFeature {
	
	public AnnotatedFeature(String nameOfChromosome, Strand strand, int coordinateStart, int coordinateEnd) {
		super(nameOfChromosome, strand, coordinateStart, coordinateEnd);
	}
	
    private Map<String, Gene> gene2transcripts = new HashMap<String, Gene>();
    	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) return false;
		if (this == obj) return true;
		if (!(obj instanceof ContiguousGenomicFeature)) return false;
		ContiguousGenomicFeature that = (ContiguousGenomicFeature)obj;
		if (!this.getType().equals(that.getType())) return false;
		if (!this.getIdOfReferenceSequence().equals(that.getIdOfReferenceSequence())) return false;
		if (this.getCoordinateOfStart() != that.getCoordinateOfStart()) return false;
		if (this.getCoordinateOfEnd() != that.getCoordinateOfEnd()) return false;
		if (this.getStrand() != that.getStrand()) return false;
		return true;
	}
    	
	@Override
	public int hashCode() {
		// TODO Auto-generated method stub
		return (this.getType()+this.getIdOfReferenceSequence()+this.getCoordinateOfStart()+this.getCoordinateOfEnd()+this.getStrand()).hashCode();
	}
	
	public boolean addTranscript(String geneId, String transcriptId) {
		if (geneId == null || transcriptId == null) return false;
		Gene gene = gene2transcripts.get(geneId);
		if (gene == null) {
			gene = new Gene();
			gene.geneId = geneId;
			gene2transcripts.put(geneId, gene);
		}
		return gene.addTransciptId(transcriptId);
	}
	
	public Set<Gene> getGenes() {
		return Collections.unmodifiableSet(new HashSet<Gene>(gene2transcripts.values()));
	}
	
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		StringBuffer buffer = new StringBuffer(super.toString());
		buffer.deleteCharAt(buffer.length()-1);
		for (Map.Entry<String, Gene> entry : gene2transcripts.entrySet()) {
			buffer.append(String.format("gene_id \"%s\";", entry.getKey()));
			for (String transcriptId : entry.getValue().getTranscriptIds()) {
				buffer.append(String.format(" transcript_id \"%s\";", transcriptId));
			}
		}
		return buffer.toString();
	}
	
	public static class Gene {
		private String geneId;
		private Set<String> transcriptIds = new HashSet<String>();
		private Gene() {}
			
		public String getId() {
			return geneId;
		}
		
		public Set<String> getTranscriptIds() {
			return Collections.unmodifiableSet(transcriptIds);
		}
		
		boolean addTransciptId(String transcriptId) {
			if (transcriptId == null) return false;
			return transcriptIds.add(transcriptId);
		}
	}
};
