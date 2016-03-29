package com.lifetechnologies.solid.wt.splice;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.ParseException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;


import com.lifetechnologies.solid.wt.ContiguousGenomicFeature;
import com.lifetechnologies.solid.wt.ContiguousGenomicFeatureComparator;
import com.lifetechnologies.solid.wt.GFFFileUtilities;
import com.lifetechnologies.solid.wt.IndexedFastaDatabase;
import com.lifetechnologies.solid.wt.InvalidIndexException;
import com.lifetechnologies.solid.wt.Strand;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.WTException;
import com.lifetechnologies.util.BaseStream;

/**
 * Extract the sequence flanking known and putative splice junctions.
 * A known splice junction is a junction between exons in a known transcript.
 * A putative splice junction is a junction between exons in a gene that is not
 * in a known transcript.
 * 
 * @author mullermw
 *
 */
public class SpliceJunctionExtractor {

	private Logger logger = Logger.getLogger(SpliceJunctionExtractor.class.toString());
	private IndexedFastaDatabase sequenceReference;
	private InputStream exonReferenceGtf;
	private PrintStream out;
	private int flankSize;
	
	/**
	 * 
	 * @param sequenceReference database of reference sequence.
	 * @param exonReferenceGff a GTF file containing the exons in the reference.
	 * @param out Where to write the extracted junctions.
	 * @param flankSize Maximum amount of flanking sequence to extract.
	 */
	public SpliceJunctionExtractor(IndexedFastaDatabase sequenceReference, InputStream exonReferenceGff, PrintStream out, int flankSize) {
		if (sequenceReference == null) throw new IllegalArgumentException("sequenceReference cannot be null");
		if (exonReferenceGff == null) throw new IllegalArgumentException("exonReferenceGff cannot be null");
		if (out == null) out = System.out;
		if (flankSize < 0) throw new IllegalArgumentException("flankSize cannot be less than 1");
		this.sequenceReference = sequenceReference;
		this.exonReferenceGtf = exonReferenceGff;
		this.out = out;
		this.flankSize = flankSize;
	}
	
	/**
	 * Extract the junctions.
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 * @throws InvalidIndexException
	 * @throws WTException
	 */
	public long extractJunctions() throws IOException, ParseException, InvalidIndexException, WTException {
		logger.finer("Parsing GTF file");
		List<AnnotatedFeature> features = GFFFileUtilities.parseUniqueExonsFromGtf(exonReferenceGtf);
		//contig -> geneid -> Gene
		Map<String, Map<String, Gene>> geneIndex = new HashMap<String, Map<String, Gene>>();
		Collection<String> referenceHeaders = sequenceReference.getListOfHeaders(false/*don't store*/, true/*truncateAfterFirstSpace*/);
		for (String contigId : referenceHeaders)
			geneIndex.put(contigId, new HashMap<String, Gene>());
		Set<String> ignoredContigIds = new HashSet<String>();
		for (AnnotatedFeature feature : features) {

			String contigId = feature.getIdOfReferenceSequence();
			if (!referenceHeaders.contains(contigId)) {
				ignoredContigIds.add(contigId);
				continue;
			}
			Map<String, Gene> geneId2Gene = geneIndex.get(contigId);
			for (AnnotatedFeature.Gene aGene : feature.getGenes()) {
				String geneId = aGene.getId();
				Gene gene = geneIndex.get(contigId).get(geneId);
				if (gene == null) {
					gene = new Gene(geneId);
					gene.setChromId(contigId);
					geneId2Gene.put(geneId, gene);
				}
				for (String transcriptId : aGene.getTranscriptIds()) {
					gene.add(feature, transcriptId);
				}
			}
		}

		if (ignoredContigIds.size() > 0)
			logger.info("Ignored exon features that don't correspond to the reference:"+ignoredContigIds);
		ignoredContigIds = null;
		
		int junctionCounter = 0;
		//Loop once per contig.
		for (String contig : referenceHeaders) {
			logger.info("Extracting Junctions on " + contig);
			List<Gene> genes = new ArrayList<Gene>(geneIndex.get(contig).values());
			
			List<FlankedJunction> junctions = new ArrayList<FlankedJunction>();
			for (Gene gene : genes) 
				junctions.addAll(gene.getJunctions(flankSize));
	
			final List<String> headerList = sequenceReference.getListOfHeaders(false, true);
			//Sort the junctions location in the sequence reference file.
			Collections.sort(junctions, new Comparator<FlankedJunction>() {
				@Override
				public int compare(FlankedJunction o1, FlankedJunction o2) {
					int idx1 = headerList.indexOf(o1.getGene().getChromId());
					int idx2 = headerList.indexOf(o2.getGene().getChromId());
					if (idx1 != idx2) return idx1 - idx2;
					return o1.compareTo(o2);
				}
			});
	
			Deque<FlankedJunction> junctionDeque = new ArrayDeque<FlankedJunction>(junctions);
			junctions = null; //To be garbage collected.
			List<FlankedJunction> junctionsExtendingOnLeft = new ArrayList<FlankedJunction>();
			final NavigableSet<FlankedJunction> halfJunctions = new TreeSet<FlankedJunction>(new Comparator<FlankedJunction>() {
				@Override
				public int compare(FlankedJunction o1, FlankedJunction o2) {
					int idx1 = headerList.indexOf(o1.getGene().getChromId());
					int idx2 = headerList.indexOf(o2.getGene().getChromId());
					if (idx1 != idx2) return idx1 - idx2;
					if (o1.getRightStart() != o2.getRightStart()) return o1.getRightStart() - o2.getRightStart();
					return o1.compareTo(o2);
				}
			});
			List<FlankedJunction> junctionsExtendingOnRight = new LinkedList<FlankedJunction>();
			//BaseStream stream = sequenceReference.nextSequenceAsStream();
			BaseStream stream = sequenceReference.getSequenceAsStream(contig);
			try {
				while (stream.hasMoreBases() && !(junctionDeque.isEmpty() && junctionsExtendingOnLeft.isEmpty() && halfJunctions.isEmpty() && junctionsExtendingOnRight.isEmpty())) {
					FlankedJunction junction = junctionDeque.peekFirst();
					Character nextBase = stream.nextBase();
					
					//Arrived at the first base of a junction?  Move it to the left extension list.
					while (junction != null && stream.getBaseCount() == junction.getLeftStart()) {
						junctionsExtendingOnLeft.add(junctionDeque.pollFirst());
						junction = junctionDeque.peekFirst();
					}
	
					//Extend the left sides of junctions.
					for (Iterator<FlankedJunction> it=junctionsExtendingOnLeft.iterator(); it.hasNext();) {
						junction = it.next();
						junction.getLeftSequence().append(nextBase);
						if (junction.getLeftEnd() == stream.getBaseCount() ) {
							it.remove();
							halfJunctions.add(junction);
						}
					}
	
					//Arrived at the first base of the right side of a junction?  Move it to the right extension list.
					while (!halfJunctions.isEmpty() && stream.getBaseCount() == halfJunctions.first().getRightStart()) {
						junctionsExtendingOnRight.add(halfJunctions.pollFirst());
					}
					
					//Extend the right sides of junctions.
					for (Iterator<FlankedJunction> it=junctionsExtendingOnRight.iterator(); it.hasNext();) {
						junction = it.next();
						junction.getRightSequence().append(nextBase);
						if (junction.getRightEnd() == stream.getBaseCount()) {
							it.remove();
							junctionCounter++;
							out.println(">junction"+junctionCounter+" "+junction.toString());
							out.println(junction.getSequence());
						}
					}
				}
			} finally {
				stream.close();
			}
			if (!(junctionDeque.isEmpty() && junctionsExtendingOnLeft.isEmpty() && halfJunctions.isEmpty() && junctionsExtendingOnRight.isEmpty())) {
				logger.info(Utilities.join(":", stream.getBaseCount(), junctionDeque.size(), junctionsExtendingOnLeft.size(), halfJunctions.size(), junctionsExtendingOnRight.size(), stream.hasMoreBases()));
				throw new WTException("Some splice junctions could not be retrieved.");
			}
			logger.info("Extracted "+junctionCounter+" junctions from "+contig);
		}
		return junctionCounter;
	}
	
	public static void main(String[] args) throws Exception  {
		String fastaFilename = args.length > 0 ? args[0] : "test/input/human_chr17_6.fa";
		String indexDir = new File("ignore").exists() ? "ignore/myindex" : "myindex";
		String gffFilename = args.length > 1 ? args[1] : "test/input/human_chr17_6.exons.newformat.gtf";
		//String fastaFilename = args.length > 0 ? args[0] : "test/ignore/human.fa";
		//String indexDir = new File("ignore").exists() ? "ignore/human.fa.idx" : "myindex";
		//String gffFilename = args.length > 1 ? args[1] : "test/input/human_chr17_6.exons.newformat.gff";
		IndexedFastaDatabase db = new IndexedFastaDatabase(new java.io.File(fastaFilename), new java.io.File(indexDir));	
		new SpliceJunctionExtractor(db, new FileInputStream(gffFilename),System.out, 50).extractJunctions();
	}
}

/**
 * Supporting class for SpliceJunctionExtractor.
 * Represents a Gene as a set of Transcript.
 * Only the coordinates are stored, no sequence.
 * @author mullermw
 *
 */
class Gene extends TreeSet<ContiguousGenomicFeature> implements Comparable<Gene> {

	public static final long serialVersionUID = 1;
	
	private String chromId;
	private int start = Integer.MAX_VALUE;
	private int end = Integer.MIN_VALUE;
	private Strand strand;
	private String id = "";
	private Logger logger = Logger.getLogger(Gene.class.toString());
	
	private Map<String, Transcript> transcripts = new HashMap<String, Transcript>();
	
	Gene(String id ) {
		super(ContiguousGenomicFeatureComparator.INSTANCE);
		this.id = id;
	}
	
	public int getStart() {
		return start;
	}
	
	public boolean add(ContiguousGenomicFeature e, String transcriptId) {
		if (e == null) return false;
		if (transcriptId == null) throw new IllegalArgumentException("transcriptId cannot be null.");
		if (Utilities.isBlank(transcriptId)) logger.warning("transcriptId is blank.");
		if (chromId == null) chromId = e.getIdOfReferenceSequence();
		if (strand == null) strand = e.getStrand();
		if (e.getCoordinateOfStart() < start) start = e.getCoordinateOfStart();
		if (e.getCoordinateOfEnd() > end) end = e.getCoordinateOfEnd();
		Transcript transcript = transcripts.get(transcriptId);
		boolean newTranscript = false;
		if (transcript == null) {
			transcript = new Transcript(transcriptId, this);
			transcripts.put(transcriptId, transcript);
			newTranscript = true;
		}
		transcript.add(e);
		this.add(e);
		return newTranscript;
	}

	public SortedSet<FlankedJunction> getJunctions(int flankSize) {
		SortedSet<FlankedJunction> junctions = new TreeSet<FlankedJunction>();
		for (Transcript transcript : transcripts.values()) {
			try {
				junctions.addAll(transcript.getJunctions(flankSize));
			} catch (IllegalArgumentException e) {
				System.err.println(this.toString());
				System.exit(1);
				logger.severe("Error processing: "+this.toString());
				throw e;
			}
		}
		List<ContiguousGenomicFeature> exons = new ArrayList<ContiguousGenomicFeature>(this);
		for (int i=0; i<exons.size(); i++ ) {
			ContiguousGenomicFeature leftExon = exons.get(i);
			for (int j=i+1; j<exons.size(); j++) {
				ContiguousGenomicFeature rightExon = exons.get(j);
				if (leftExon.getCoordinateOfEnd() >= rightExon.getCoordinateOfStart()) continue;
				FlankedJunction junction = new FlankedJunction(this);
				junction.setKnown(false);
				junction.setLeftStart(leftExon.getCoordinateOfEnd() - flankSize + 1 > leftExon.getCoordinateOfStart() ? leftExon.getCoordinateOfEnd() - flankSize + 1 : leftExon.getCoordinateOfStart());
				junction.setLeftEnd(leftExon.getCoordinateOfEnd());
				junction.setRightStart(rightExon.getCoordinateOfStart());
				junction.setRightEnd(rightExon.getCoordinateOfStart() + flankSize - 1 < rightExon.getCoordinateOfEnd() ? rightExon.getCoordinateOfStart() + flankSize - 1: rightExon.getCoordinateOfEnd());
				junctions.add(junction);
			}
		}
		return junctions;
	}
	
	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}
	
	public String getChromId() {
		return chromId;
	}

	public void setChromId(String chromId) {
		this.chromId = chromId;
	}

	public Strand getStrand() {
		return strand;
	}

	public void setStrand(Strand strand) {
		this.strand = strand;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	@Override
	public boolean equals(Object o) {
		if (o == null) return false;
		if (this == o ) return true;
		Gene that = (Gene)o;
		if (this.compareTo(that) == 0) return true;
		return false;
	} 
	
	@Override
	public int compareTo(Gene that) {
		if (that == null) return -1;
		if (this == that) return 0;
		if (!this.getChromId().equals(that.getChromId()))
			return this.getChromId().compareTo(that.getChromId());
		if (this.getStart() != that.getStart()) 
			return this.getStart() - that.getStart();
		if (this.getEnd() != that.getEnd()) 
			return this.getEnd() - that.getEnd();
		return this.getStrand().compareTo(that.getStrand());
	}
	
	public String toString() {
		StringBuffer val = new StringBuffer();
		val.append("Gene:"+id+"["+chromId+","+start+"-"+end+"]\n");
		for (ContiguousGenomicFeature feature : this) {
			val.append("   "+feature+"\n");
		}
		for (Transcript transcript : transcripts.values()) {
			val.append("   "+transcript+"\n");
		}
		return val.toString();
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
}

/**
 * Supporting Class for SpliceJunctionExtractor
 * Represents a Transcript as a set of exons.
 * Only stores coordinates, not sequence.
 * @author mullermw
 *
 */
class Transcript extends TreeSet<ContiguousGenomicFeature> implements Comparable<Transcript> {
	
	public static final long serialVersionUID = 1;
	
	private int start = Integer.MAX_VALUE;
	private int end = Integer.MIN_VALUE;
	private String id = "";
	private Gene gene;
	
	Transcript(String id, Gene gene) {
		super(ContiguousGenomicFeatureComparator.INSTANCE);
		this.gene = gene;
		this.id = id;
	}
	
	@Override
	public boolean add(ContiguousGenomicFeature e) {
		if (e == null) return false;
		if (e.getCoordinateOfStart() < start) start = e.getCoordinateOfStart();
		if (e.getCoordinateOfEnd() > end) end = e.getCoordinateOfEnd();
		return super.add(e);
	}
	
	public List<FlankedJunction> getJunctions(int flankSize) {
		List<FlankedJunction> junctions = new ArrayList<FlankedJunction>();
		ContiguousGenomicFeature precedingExon = null;
		for (ContiguousGenomicFeature exon : this) {
			if (precedingExon != null) {
				FlankedJunction junction = new FlankedJunction(gene);
				junction.setKnown(true);
				junction.setLeftStart(precedingExon.getCoordinateOfEnd() - flankSize + 1> precedingExon.getCoordinateOfStart() ? precedingExon.getCoordinateOfEnd() - flankSize + 1 : precedingExon.getCoordinateOfStart());
				junction.setLeftEnd(precedingExon.getCoordinateOfEnd());
				junction.setRightStart(exon.getCoordinateOfStart());
				junction.setRightEnd(exon.getCoordinateOfStart() + flankSize - 1 < exon.getCoordinateOfEnd() ? exon.getCoordinateOfStart() + flankSize - 1 : exon.getCoordinateOfEnd());
				junctions.add(junction);
			}
			precedingExon = exon;
		}
		return junctions;
	}
	
	@Override
	public int compareTo(Transcript that) {
		if (that == null) return -1;
		if (this == that) return 0;
		if (this.getStart() != that.getStart()) 
			return this.getStart() - that.getStart();
		return this.getEnd() - that.getEnd();
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}
	public String toString() {
		StringBuffer val = new StringBuffer();
		val.append("Transcript:"+id+"["+start+"-"+end+"]\n");
		for (ContiguousGenomicFeature feature : this) {
			val.append("      "+feature+"\n");
		}
		return val.toString();
	}
}