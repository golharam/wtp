package com.lifetechnologies.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.ResourceBundle;
import java.util.Set;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.Strand;
import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.config.WTConfigUtils;
import com.lifetechnologies.solid.wt.tile.SimpleTile;
import com.lifetechnologies.solid.wt.tile.Tile;
import com.lifetechnologies.solid.wt.tile.TileUtils;

/**
 * 
 * Reads a refGene formatted stream and writes a GTF formatted one from it.
 * 
 * refGene.txt is a file one can download from the UCSC Genome browser.
 * http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz
 * 
 * This file is a tab delimited dump from a relational database table:
 * 
 * CREATE TABLE `refGene` (
 * `bin` smallint(5) unsigned NOT NULL default '0',
 * `name` varchar(255) NOT NULL default '',
 * `chrom` varchar(255) NOT NULL default '',
 * `strand` char(1) NOT NULL default '',
 * `txStart` int(10) unsigned NOT NULL default '0',
 * `txEnd` int(10) unsigned NOT NULL default '0',
 * `cdsStart` int(10) unsigned NOT NULL default '0',
 * `cdsEnd` int(10) unsigned NOT NULL default '0',
 * `exonCount` int(10) unsigned NOT NULL default '0',
 * `exonStarts` longblob NOT NULL,
 * `exonEnds` longblob NOT NULL,
 * `id` int(10) unsigned NOT NULL default '0',
 * `name2` varchar(255) NOT NULL default '',
 * `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL default 'none',
 * `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL default 'none',
 * `exonFrames` longblob NOT NULL,
 * KEY `chrom` (`chrom`(7),`bin`),
 * KEY `name2` (`name2`(10)),
 * KEY `name` (`name`(12)),
 * KEY `chrom_2` (`chrom`(7),`txStart`),
 * KEY `chrom_3` (`chrom`(7),`txEnd`)
 * ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
 *
 * GTF Specification here: http://mblab.wustl.edu/GTF22.html
 * 
 * While a GTF export of refGene is available from UCSC, it does not contain
 * the geneID information we need.
 * 
 * When a transcript_id (refGene.name) is repeated within a contig (refGene.chrom)
 * the transcript_id will have a _dup[0-9]+ subscript added to indicate that the
 * mapping is ambiguous.
 * 
 * When a gene_id (refGene.name2) is repeated in non-overlapping locations 
 * within a contig (refGene.chrom) the gene_id will have a _dup[0-9]+ subscript
 * added to indicate that the mapping is ambigous.
 * 
 * @author mullermw
 *
 */
public class RefGene2Gtf {

	private InputStream refGene;
	private PrintStream gtf;
	private Properties header = new Properties();
	private static Logger logger = Logger.getLogger(RefGene2Gtf.class.toString());
	public static final String SOURCE_FILE_HEADER_KEY = "source-file";
	
	/**
	 * 
	 * @param refGene the refGene content to read.
	 * @param gtf destination stream for the gtf.
	 * @param header
	 */
	public RefGene2Gtf(InputStream refGene, PrintStream gtf, Properties header) {
		if (refGene == null) throw new IllegalArgumentException("refGene cannot be null.");
		if (gtf == null) this.gtf = System.out;
		if (header == null) header = new Properties();
		this.refGene = refGene;
		this.gtf = gtf;
		this.header = header;
	}
	
	/**
	 * Execute the conversion.
	 * @throws IOException
	 */
	public void run() throws IOException {
		gtf.println("## gff-version 2");
		gtf.println("## gtf");
		ResourceBundle resourceBundle = WTConfigUtils.getApplicationProperties();
		gtf.println("## source-version refgene2gff.sh " + resourceBundle.getString("version"));
		for (Map.Entry<Object, Object> entry: header.entrySet())
			gtf.printf("## %s %s\n", entry.getKey(), entry.getValue());
		gtf.println("## date " + new SimpleDateFormat("yyyy-MM-dd").format(new Date()));
		gtf.println("## This file is a transformation of the refGene.txt file from the ");
		gtf.println("## UCSC genome browser FTP site.");
		gtf.println("## Example: ");
		gtf.println("## http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz ");
		gtf.println("##");
		gtf.println("## This file is similar to the gtf files available through the UCSC genome");
		gtf.println("## browser, with a few differnces.");
		gtf.println("## The gene_id attribute contains the HUGO name of the associated gene.");
		gtf.println("## Note that when a gene_id is repeated at multiple loci in a contig/chromosome,");
		gtf.println("## a subscript will be added in the following form: GENE_dup1, GENE_dup2, etc.");
		gtf.println("## Note that when a transcript_id is repeated at multiple loci in a contig/");
		gtf.println("## chromosome, a subscript will be added in the following form:");
		gtf.println("## TRANSCRIPT_dup1, TRANSCRIPT_dup2, etc.");
		gtf.println("##");
		gtf.println("##");
		BufferedReader reader = null;
		
		//As we read in the transcripts, we aggregate them by geneId in 'geneId2Transcripts'.
		Map<String, Set<Transcript>> geneId2Transcripts = new HashMap<String, Set<Transcript>>();
		
		try {
			reader = new BufferedReader(new InputStreamReader(refGene));
			int linecount = 0;
			for (String line = reader.readLine(); line != null; line=reader.readLine()) {
				linecount++;
				line = line.trim();
				if (line.startsWith("#")) continue;
				String[] fields = line.split("\t");
				if (fields.length != 16) {
					logger.warning("Line "+linecount+" does not match refGene format: "+line);
					continue;
				}
				//Parse a row and populate a Transcript object.
				String chr = fields[2];
				String strand = fields[3];
				int geneStart = Integer.parseInt(fields[4]) + 1;
				int geneEnd = Integer.parseInt(fields[5]);
				Transcript transcript = new Transcript(chr, Strand.toStrand(strand.charAt(0)), geneStart, geneEnd);
				for (String s : fields[9].split("\\s*,\\s*"))
					transcript.exonStarts.add(Integer.parseInt(s) + 1);
				for (String s : fields[10].split("\\s*,\\s*"))
					transcript.exonEnds.add(Integer.parseInt(s));
				transcript.geneId = fields[12];
				transcript.transcriptId = fields[1];
				transcript.cdsStart = Integer.parseInt(fields[6]) + 1;
				transcript.cdsEnd = Integer.parseInt(fields[7]);
				transcript.cdsStartStat = fields[13];
				transcript.cdsEndStat = fields[14];
				for (String s : fields[15].split("\\s*,\\s*"))
					transcript.phases.add(Integer.parseInt(s));
				//Index the transcript by geneId
				Set<Transcript> transcripts = geneId2Transcripts.get(transcript.geneId);
				if (transcripts == null ) {
					transcripts = new HashSet<Transcript>();
					geneId2Transcripts.put(transcript.geneId, transcripts);
				}
				transcripts.add(transcript);
			}
		} finally {
			if (reader != null) reader.close();
		}

		//Deal with every GeneId independently.
		//As transcripts are processed, they are removed from geneId2Transcripts
		//and added to allTranscripts for downstream processing.
		List<Transcript> allTranscripts = new ArrayList<Transcript>();
		for (Iterator<Map.Entry<String, Set<Transcript>>> it = geneId2Transcripts.entrySet().iterator(); it.hasNext();) {
			Map.Entry<String, Set<Transcript>> entry = it.next();
			it.remove();
			
			//Only one transcript.
			if (entry.getValue().size() < 2) {
				allTranscripts.addAll(entry.getValue());
				continue;
			}
			
			List<Transcript> transcripts = new ArrayList<Transcript>(entry.getValue());
			Collections.sort(transcripts);
			//Count the occurences of a transcript_id per contig.
			//Contig->Transcript->Counter
			Map<String, Map<String, Counter>> transcriptCountingMatrix = new HashMap<String, Map<String, Counter>>();
			for (Transcript transcript : transcripts) {
				Map<String, Counter> transcriptId2Counter = transcriptCountingMatrix.get(transcript.getContigId());
				if (transcriptId2Counter == null) {
					transcriptId2Counter = new HashMap<String, Counter>();
					transcriptCountingMatrix.put(transcript.getContigId(), transcriptId2Counter);
				}
				if (transcriptId2Counter.containsKey(transcript.transcriptId))
					transcriptId2Counter.get(transcript.transcriptId).increment();
				else
					transcriptId2Counter.put(transcript.transcriptId, new Counter(1));
			}

			//Assign dup## subscripts to transcript_ids that appear multiple times in a contig.
			for (Map.Entry<String, Map<String, Counter>> contigCounts : transcriptCountingMatrix.entrySet()) {
				String contigId = contigCounts.getKey();
				for (String transcriptId : contigCounts.getValue().keySet()) {
					if (contigCounts.getValue().get(transcriptId).value < 2) continue;
					int counter = 0;
					for (Transcript transcript : transcripts)
						if (transcript.getContigId().equals(contigId) &&
							transcript.transcriptId.equals(transcriptId))
							transcript.transcriptId += "_dup" + (++counter);
				}
			}

			// Cluster the transcripts to determine the number of non-overlapping 'genes'
			Set<Set<Transcript>> clusters = TileUtils.clusterTiles(transcripts);
			if (clusters.size() < 2) {
				allTranscripts.addAll(transcripts);
				continue;
			}

			//sort the clusters by lowest genomic location.
			List<Set<Transcript>> clusterList = new ArrayList<Set<Transcript>>(clusters);
			Collections.sort(clusterList, new Comparator<Set<Transcript>>() {
				@Override
				public int compare(Set<Transcript> o1, Set<Transcript> o2) {
					List<Transcript> list1 = new ArrayList<Transcript>(o1);
					List<Transcript> list2 = new ArrayList<Transcript>(o2);
					Collections.sort(list1);
					Collections.sort(list2);
					return list1.get(0).compareTo(list2.get(0));
				}
			});

			//Count the non-overlapping occurences of a gene id per contig.
			Map<String, Counter> contig2geneIdCount = new HashMap<String, Counter>();
			boolean overlapFound = false;
			for (Set<Transcript> cluster : clusterList) {
				String contigId = cluster.iterator().next().getContigId();
				Counter countForThisContig = contig2geneIdCount.get(contigId);
				if (countForThisContig == null) {
					countForThisContig = new Counter(1);
					contig2geneIdCount.put(contigId, countForThisContig);
				} else {
					countForThisContig.increment();
					overlapFound = true;
				}
			}
			
			//No overlaps.
			if (!overlapFound) {
				allTranscripts.addAll(transcripts);
				continue;
			}
			
			//Assign _dup## subscripts to genes that repeat.
			Map<String, Counter> contig2geneIdDupNum = new HashMap<String, Counter>();
			for (Set<Transcript> cluster : clusterList) {
				String contigId = cluster.iterator().next().getContigId();
				int countForThisContig = contig2geneIdCount.get(contigId).value;
				if (countForThisContig > 1) {
					Counter counter = contig2geneIdDupNum.get(contigId);
					if (counter == null) {
						counter = new Counter(0);
						contig2geneIdDupNum.put(contigId, counter);
					}
					counter.increment();
					for (Transcript transcript : cluster) {
						transcript.geneId += "_dup"+counter.value;
					}
				}
				allTranscripts.addAll(cluster);
			}
		}
		
		// Sort the transcripts by Genomic location.
		Collections.sort(allTranscripts);
		
		//Write the Transcripts to the GTF file.
		for (Transcript transcript : allTranscripts) {
			for (int i=0; i<transcript.exonStarts.size(); i++) {
				Integer exonStart = transcript.exonStarts.get(i);
				Integer exonEnd = transcript.exonEnds.get(i);
				Strand strand = transcript.getStrand();
				Tile exon = SimpleTile.newTile(transcript.getContigId(), transcript.getStrand(), exonStart, exonEnd);
				Integer frame = transcript.phases.get(i);
				//Note that frame/phase has different meaning in refGene than in GTF.
				//In refGene, frame is the position of the first base in the
				//codon.  In GTF, phase is the first position in the CDS that
				//is the first base in a codon.  So we make the following
				//conversion.
				if (frame == 1) frame = 2;
				else if (frame == 2) frame = 1;
				String phase = frame > -1 ? frame.toString() : ".";
				String attributes = "gene_id \""+transcript.geneId+"\"; transcript_id \""+transcript.transcriptId+"\";";
				int cdsStart = exonStart;
				int cdsEnd = exonEnd;
				boolean isFirstCodingExon = false;
				
				//Conditionally write the start_codon feature.
				//Following the conventions of GTF exporting at UCSC.
				if (strand == Strand.POSITIVE) {
					isFirstCodingExon = !transcript.cdsStartStat.equals("unk") && exon.contains(transcript.cdsStart);
				} else {
					isFirstCodingExon = !transcript.cdsEndStat.equals("unk") && exon.contains(transcript.cdsEnd);
				}
				if (isFirstCodingExon) {
					int startCodonStart = strand == Strand.POSITIVE ? transcript.cdsStart : transcript.cdsEnd - 2;
					gtf.println(Utilities.join("\t", 
							transcript.getContigId(),
							"refGene",
							"start_codon",
							startCodonStart,
							startCodonStart + 2,
							"0.000000",
							transcript.getStrand().toChar(),
							".",
							attributes
						));
				}

				//Conditionally write the CDS feature.
				if (transcript.cdsStartStat.equals("cmpl") && exon.contains(transcript.cdsStart))
					cdsStart = transcript.cdsStart;
				if (transcript.cdsEndStat.equals("cmpl") && exon.contains(transcript.cdsEnd))
					cdsEnd = transcript.cdsEnd;

				if (!phase.equals("."))
					gtf.println(Utilities.join("\t", 
							transcript.getContigId(),
							"refGene",
							"CDS",
							cdsStart,
							cdsEnd,
							"0.000000",
							transcript.getStrand().toChar(),
							phase,
							attributes
					));
				
				//Write the exon feature.
				gtf.println(Utilities.join("\t", 
						transcript.getContigId(), 
						"refGene", 
						"exon", 
						exonStart, 
						exonEnd, 
						"0.000000",
						transcript.getStrand().toChar(),
						".",
						attributes
						));
			}
		}
	}
	
	/**
	 * Represents a Transcript, contains fields corresponding to a row in the refGene.txt file.
	 * @author mullermw
	 *
	 */
	class Transcript extends SimpleTile {
		String geneId;
		String transcriptId;
		String cdsStartStat;
		String cdsEndStat;
		int cdsStart;
		int cdsEnd;
		//int geneIdInstance = 0;
		//int transcriptIdInstance = 0;
		List<Integer> exonStarts = new ArrayList<Integer>();
		List<Integer> exonEnds = new ArrayList<Integer>();
		List<Integer> phases = new ArrayList<Integer>();
		
		public Transcript(String contigId, Strand strand, int start, int end) {
			super(contigId, strand, start, end);
		}
		
		@Override
		protected Object clone() throws CloneNotSupportedException {
			Transcript clone = new Transcript(this.getContigId(), this.getStrand(), this.getStart(), this.getEnd());
			clone.geneId = this.geneId;
			clone.transcriptId = this.transcriptId;
			clone.exonStarts = new ArrayList<Integer>(exonStarts);
			clone.exonEnds = new ArrayList<Integer>(exonEnds);
			clone.phases = new ArrayList<Integer>(phases);
			return clone;
		}
		
		@Override
		public String toString() {
			// TODO Auto-generated method stub
			return super.toString()+this.geneId+":"+this.transcriptId;
		}
	}
	
	/**
	 * Simple substitute for an int pointer.
	 * @author mullermw
	 *
	 */
	class Counter {
		int value = 0;
		
		public Counter(int initialValue) {
			value = initialValue;
		}
		
		public int  increment() {
			return ++value;
		}	
	}
	
}


