package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import com.lifetechnologies.util.SeekableIterator;
import com.lifetechnologies.util.SeekablePolyIterator;

/**
 * User: tuchbb
 * Date: Nov 19, 2008
 * Time: 11:38:30 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 * 
 * This is the fundamental exon finding code.  It will predict Transcribed Regions(TR) using the reads
 * from specified region(s) of the reference.
 * 
 */
public class ExonFinder {

	public static final String AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON = "avg-frac-ptr-overlapping-atr";
	public static final String AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON = "avg-frac-atr-overlapping-ptr";
	public static final String COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS = "number-of-atrs-matching-ptr";
	public static final String COUNT_OF_ALL_PREDICTED_EXONS_OVERLAPPING_AN_ANNOTATED_EXON = "number-of-ptrs-matching-atr";
	public static final String COUNT_OF_ALL_ANNOTATED_EXONS = "atr-count";
	public static final String COUNT_OF_ALL_PREDICTED_EXONS = "ptr-count";
	public static final String MIN_FRACTION_OF_AVERAGE_COVERAGE_BEFORE_TRIMMING_POSITION_FROM_PREDICTED_EXON = "trimming-fraction";
	public static final String MIN_ALIGNMENT_SCORE = "min-score";
	public static final String MIN_OVERLAP_FOR_COUNTING_PREDICTED_EXON_AS_TRUE_POSITIVE = "min-overlap";
	public static final String MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON = "min-window-coverage";
	public static final String WIDTH_OF_MOVING_WINDOW = "window-size";
	private static final Logger logger = Logger.getLogger(ExonFinder.class.toString());
	private static boolean showProgress = false;
	public static final String INDEX_DIR_SYS_PROPERTY = "com.lifetechnologies.solid.wt.FastaIndexDir";
	public static final String SHOW_PROGRESS_SYS_PROPERTY = "com.lifetechnologies.solid.wt.showProgress";
	private static final NumberFormat COVERAGE_FORMAT = new DecimalFormat("0.0");
	private static final NumberFormat OVERLAP_FORMAT = new DecimalFormat("0.00");
	
	
	public static void main(String args[]) throws Exception {

		if (args.length < 9 || args.length > 10) {
			System.err.println(
					"Usage:\n" +
					"ExonFinder "+WIDTH_OF_MOVING_WINDOW+"\n" +
					"           "+MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON+"\n"+
					"           "+MIN_OVERLAP_FOR_COUNTING_PREDICTED_EXON_AS_TRUE_POSITIVE+"\n"+
					"           "+MIN_ALIGNMENT_SCORE+"\n"+
					"           "+MIN_FRACTION_OF_AVERAGE_COVERAGE_BEFORE_TRIMMING_POSITION_FROM_PREDICTED_EXON+"\n"+
					"           fileListOfSortedMAXFiles\n"+
					"           fastaDatabaseOfReference\n"+
					"           fileExonAnnotationGFF\n"+
					"           folderForOutput\n"+
					"           [comma separated list of Genomic locations of the form <refidx>:<start>-<end>] default: all sequences in reference\n"+
					"\n"+
					" use -D"+INDEX_DIR_SYS_PROPERTY+"=PATH to specify the path to the fasta file index.\n"+
					" use -D"+SHOW_PROGRESS_SYS_PROPERTY+"=false to supress the progress message.\n"
			);
			System.exit(1);
		}

		int widthOfMovingWindow = Integer.parseInt(args[0]);
		double minCoverageToCallAsPartOfExon = Double.parseDouble(args[1]);
		double minOverlapForCountingPredictedExonAsTruePositive = Double.parseDouble(args[2]);
		int minAlignmentScore = Integer.parseInt(args[3]);
		double minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon = Double.parseDouble(args[4]);
		File fileListOfSortedMAXFiles = new File(args[5]);
		File referenceFastaFile = new File(args[6]);
		File fileExonAnnotationGFF = Utilities.isBlank(args[7]) ? null : new File(args[7]);
		File folderForOutput = new File(args[8]);
		showProgress = Boolean.parseBoolean(System.getProperty(SHOW_PROGRESS_SYS_PROPERTY, "true"));
		
		//Client can change the location of the index dir with a System property.
		//Default is in the output folder.
		String indexDirString = System.getProperty(INDEX_DIR_SYS_PROPERTY);
		File indexDir = indexDirString == null ? null : new File(indexDirString);
		if (indexDir == null) {
			indexDir = folderForOutput;
		} else {
			if (indexDir.isDirectory() == false) throw new Exception("Index Dir: "+indexDir+" is not a directory");
			if (indexDir.canRead() == false) throw new Exception("Cannot read Index Dir: "+indexDir);
		}
		
		File referenceIndexFile = new File(indexDir, referenceFastaFile.getName().concat(".idx"));
		logger.info("new IndexedFastaDatabase()");
		
		//Do not overwrite index file.
		FastaDatabase fastaDatabaseOfReference = new IndexedFastaDatabase(referenceFastaFile, referenceIndexFile, !referenceIndexFile.exists());

		List<String> listOfHeaders = fastaDatabaseOfReference.getListOfHeaders(true, false);
		
		//Map reference ids back to indices
		final Map<String, Integer> referenceIdToIndex = new HashMap<String, Integer>();
		for (int i=0; i<listOfHeaders.size(); i++) {
			String id = listOfHeaders.get(i).replaceAll("\\s.*","");
			referenceIdToIndex.put(id, i);
		}

		//Slightly modified version of ContiguousGenomicFeatureComparator
		//Sorts by reference Index rather than reference id.
		Comparator<ContiguousGenomicFeature> analysisRegionComparator = new Comparator<ContiguousGenomicFeature>() {
			public int compare(ContiguousGenomicFeature contiguousGenomicFeatureA, ContiguousGenomicFeature contiguousGenomicFeatureB) {
				int refIdxA = referenceIdToIndex.get(contiguousGenomicFeatureA.getIdOfReferenceSequence());
				int refIdxB = referenceIdToIndex.get(contiguousGenomicFeatureB.getIdOfReferenceSequence());
				if (!contiguousGenomicFeatureA.equals(contiguousGenomicFeatureB)) {
					if (refIdxA == refIdxB) {
						if (contiguousGenomicFeatureA.getCoordinateOfStart() == contiguousGenomicFeatureB.getCoordinateOfStart()) {
							if (contiguousGenomicFeatureA.getCoordinateOfEnd() == contiguousGenomicFeatureB.getCoordinateOfEnd()) {
								if (contiguousGenomicFeatureA.getStrand() == contiguousGenomicFeatureA.getStrand()) {
									if (contiguousGenomicFeatureA.getType() == null || contiguousGenomicFeatureB.getType() == null
											|| contiguousGenomicFeatureA.getType().equalsIgnoreCase(contiguousGenomicFeatureB.getType())) {
										return 0;
									} else
										return contiguousGenomicFeatureA.getType().compareTo(contiguousGenomicFeatureB.getType());
								} else if (contiguousGenomicFeatureA.getStrand() == Strand.POSITIVE )
									return -1;
								else
									return 1;
							} else if (contiguousGenomicFeatureA.getCoordinateOfEnd() < contiguousGenomicFeatureB.getCoordinateOfEnd())
								return -1;
							else
								return 1;
						} else if (contiguousGenomicFeatureA.getCoordinateOfStart() < contiguousGenomicFeatureB.getCoordinateOfStart())
							return -1;
						else
							return 1;
					} else
						return refIdxA - refIdxB;
				}
				return 0;
			}
		};

		SortedSet<ContiguousGenomicFeature> regionsToAnalyze = new TreeSet<ContiguousGenomicFeature>(analysisRegionComparator);
		if (args.length > 9)
			regionsToAnalyze.addAll(ContiguousGenomicFeature.parseFeatures(args[9]));

		SortedMap<String, Long> mapIdToSequenceLength = new TreeMap<String, Long>();
		for (Map.Entry<String, Long> entry : fastaDatabaseOfReference.getSequenceLengths().entrySet()) {
			String header = entry.getKey();
			String id = header.replaceAll("\\s.*", "");
			mapIdToSequenceLength.put(id, entry.getValue());
		}
	
		//Validate genomic intervals.
		for (ContiguousGenomicFeature feature : regionsToAnalyze) {
			int sequenceLength = mapIdToSequenceLength.get(feature.getIdOfReferenceSequence()).intValue();
			if (feature.getCoordinateOfStart() < 1) feature.setCoordinateOfStart(0);
			if (feature.getCoordinateOfEnd() < 1) feature.setCoordinateOfEnd(sequenceLength+1);
			if (feature.getCoordinateOfStart() > sequenceLength - 1)
				throw new Exception("Invalid genomic interval:" + feature.toBriefString());
			if (feature.getCoordinateOfEnd() > sequenceLength - 1) {
				logger.warning("truncating genomic region: "+ feature.toBriefString());
				feature.setCoordinateOfEnd(sequenceLength - 1);
			}
			
			if (feature.getCoordinateOfStart() > feature.getCoordinateOfEnd()) {
				int tmp = feature.getCoordinateOfStart();
				feature.setCoordinateOfStart(feature.getCoordinateOfEnd());
				feature.setCoordinateOfEnd(tmp);
			}
		}
		

		//No regions specified, do everything.
		if (regionsToAnalyze.isEmpty()) {
			Set<String> ids = new HashSet<String>();
			for (String header : listOfHeaders) {
				String id = header.replaceAll("\\s.*", "");
				if (ids.contains(id)) {
					logger.warning("Ignoring duplicate sequence id: '"+id+"'.");
				} else {
					regionsToAnalyze.add(new ContiguousGenomicFeature(id, Strand.EITHER, 0, mapIdToSequenceLength.get(id).intValue() - 1));
					ids.add(id);
				} 
			}
		}
		
		System.err.println("filesMAXSorted:");
		HashSet<String> setOfMaxFilesPaths = TextFileUtilities.loadSetFromFile(fileListOfSortedMAXFiles, "\t", 0);
		//File filesMAXSorted[] = new File[setOfMaxFilesPaths.size()];
		Set<SeekableIterator<ExtendedReadMapping>> iterators = new HashSet<SeekableIterator<ExtendedReadMapping>>();

		for (Iterator<String> iteratorOverMAXFilePaths = setOfMaxFilesPaths.iterator();  iteratorOverMAXFilePaths.hasNext(); ) {
			//filesMAXSorted[indexOfFile] = new File(iteratorOverMAXFilePaths.next());
			File maxFile = new File(iteratorOverMAXFilePaths.next());
			if (!maxFile.exists()) throw new FileNotFoundException("Can't find file: " + maxFile);
			if (!maxFile.canRead()) throw new FileNotFoundException("Can't read file: " + maxFile);
			iterators.add(new SortedMaxFileIterator(maxFile));
			System.err.println(maxFile);
		}
		System.err.println();
		//SortedMAXFilesBufferedReader readerSortedMAXFiles = new SortedMAXFilesBufferedReader(filesMAXSorted, new ExtendedReadMappingByPositionComparator(), false);
		SeekablePolyIterator<ExtendedReadMapping> ermPolyIterator = new SeekablePolyIterator<ExtendedReadMapping>(ExtendedReadMappingByPositionComparator.INSTANCE,iterators);

		findExons(widthOfMovingWindow, minCoverageToCallAsPartOfExon, minOverlapForCountingPredictedExonAsTruePositive, minAlignmentScore,
				minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon,
				referenceIdToIndex,
				fileExonAnnotationGFF, ermPolyIterator, folderForOutput, regionsToAnalyze);
	}

	private static void findExons(int widthOfMovingWindow,
			double minCoverageToCallAsPartOfExon,
			double minOverlapForCountingPredictedExonAsTruePositive,
			int minAlignmentScore,
			double minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon,
			Map<String, Integer> referenceIdToIndex,
			File fileExonAnnotationGFF,
			SeekablePolyIterator<ExtendedReadMapping> ermPolyIterator,
			File folderForOutput,
			SortedSet<ContiguousGenomicFeature> regionsToAnalyze) throws Exception {
		logger.info("findExons()");

		logger.info("parsing features");
		ContiguousGenomicFeatureList listOfContiguousGenomicFeatures = null;
		if (fileExonAnnotationGFF == null) {
			listOfContiguousGenomicFeatures = new ContiguousGenomicFeatureList(1000);
		} else {
			listOfContiguousGenomicFeatures = new ContiguousGenomicFeatureList(fileExonAnnotationGFF, 1000, 0, true, true, false, true);
		}
			
		//Remove any features not contained in regionsToAnalyze and any non-exon features.
		for (Iterator<ContiguousGenomicFeature> it=listOfContiguousGenomicFeatures.getSetOfAllFeatures().iterator(); it.hasNext();) {
        	ContiguousGenomicFeature feature = it.next();
        	boolean analyze = false;
        	for (ContiguousGenomicFeature analysisRegion : regionsToAnalyze) {
        		if (analysisRegion.getIdOfReferenceSequence().equals(feature.getIdOfReferenceSequence())) {
        			if (analysisRegion.getCoordinateOfStart() == -1) analyze = true;
        			if (analysisRegion.getCoordinateOfStart() <= feature.getCoordinateOfStart()) {
        				if (analysisRegion.getCoordinateOfEnd() == -1) analyze = true;
        				if (analysisRegion.getCoordinateOfEnd() >= feature.getCoordinateOfEnd()) analyze = true;
        			}
        		} 
        	}
        	if (feature.getType().trim().toLowerCase().equals("exon") == false) analyze = false;
        	if (!analyze) it.remove();
        }
		
		//Remove any redundant exons.
		Set<FeatureKey> featureSet = new HashSet<FeatureKey>();
		for (Iterator<ContiguousGenomicFeature> it=listOfContiguousGenomicFeatures.getSetOfAllFeatures().iterator(); it.hasNext();) {
			ContiguousGenomicFeature feature = it.next();
			FeatureKey key = new FeatureKey(feature);
			if (featureSet.contains(key)) {
				it.remove();
			} else {
				featureSet.add(key);
			}
		}
		featureSet.clear();
		featureSet = null;
		
		logger.info(listOfContiguousGenomicFeatures.getSetOfAllFeatures().size()+" annotated features in genomic regions.");

		if (folderForOutput.exists()) {
			if (folderForOutput.isDirectory() == false) throw new Exception(folderForOutput + " is not a directory.");
		} else {
			if (folderForOutput.mkdir() != true) throw new Exception(folderForOutput + "/ cannot be created.");
		}

		BufferedWriter writerOutputPlusStrandGFFFile = new BufferedWriter(new FileWriter(folderForOutput.getPath() + "/exonFinder." + widthOfMovingWindow + "_" + minCoverageToCallAsPartOfExon + "_" + minOverlapForCountingPredictedExonAsTruePositive + "_" + minAlignmentScore + "_" + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon + ".plus.gff"));
		BufferedWriter writerOutputMinusStrandGFFFile = new BufferedWriter(new FileWriter(folderForOutput.getPath() + "/exonFinder." + widthOfMovingWindow + "_" + minCoverageToCallAsPartOfExon + "_" + minOverlapForCountingPredictedExonAsTruePositive + "_" + minAlignmentScore + "_" + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon + ".minus.gff"));
		try {
			writerOutputPlusStrandGFFFile.write("track name=NTRFinderPlus description=\"NTRFinder (+) "+WIDTH_OF_MOVING_WINDOW+"_" + widthOfMovingWindow + " "+MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON+"_" + minCoverageToCallAsPartOfExon + " "+MIN_ALIGNMENT_SCORE+"_" + minAlignmentScore + " "+MIN_FRACTION_OF_AVERAGE_COVERAGE_BEFORE_TRIMMING_POSITION_FROM_PREDICTED_EXON+"_" + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon + "\" useScore=1");
			writerOutputPlusStrandGFFFile.newLine();
			writerOutputMinusStrandGFFFile.write("track name=NTRFinderMinus description=\"NTRFinder (-) "+WIDTH_OF_MOVING_WINDOW+"_" + widthOfMovingWindow + " "+MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON+"_" + minCoverageToCallAsPartOfExon + " "+MIN_ALIGNMENT_SCORE+"_" + minAlignmentScore + " "+MIN_FRACTION_OF_AVERAGE_COVERAGE_BEFORE_TRIMMING_POSITION_FROM_PREDICTED_EXON+"_" + minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon + "\" useScore=1");
			writerOutputMinusStrandGFFFile.newLine();
			int countOfAllPredictedExons = 0;
			int countOfAllPredictedExonsOverlappingKnownExonOnPlusStrand = 0;
			int countOfAllPredictedExonsOverlappingKnownExonOnMinusStrand = 0;
			int countOfAllAlignmentsProcessedWithMinScore = 0;

			// exon key is refid_strand_start_end
			HashMap<String, Double> mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective = new HashMap<String, Double>();
			double sumOfFractionsOfPredictedExonOverlappingAnnotatedExon = 0.0;

			Counter readCounter = new Counter();
			long startTime = System.currentTimeMillis();
			for (ContiguousGenomicFeature analysisRegion : regionsToAnalyze) {
				String idOfCurrentReferenceSequence = analysisRegion.getIdOfReferenceSequence();
				System.err.println("Processing region: " + idOfCurrentReferenceSequence+":"+analysisRegion.getCoordinateOfStart()+"-"+analysisRegion.getCoordinateOfEnd());
				if (readCounter.getCount() > 0)
					logger.info(readCounter.getCount() / ((System.currentTimeMillis() - startTime)/1000) + " reads/s");
				int indexOfCurrentReferenceSequence = referenceIdToIndex.get(idOfCurrentReferenceSequence);
				ExtendedReadMapping prototype = new ExtendedReadMapping(null, idOfCurrentReferenceSequence, indexOfCurrentReferenceSequence+1, analysisRegion.getCoordinateOfStart(), 0,0,0,0);
				ermPolyIterator.seekTo(prototype);
				long coverageSumInWindowOnPlusStrand = 0;
				long coverageSumInWindowOnMinusStrand = 0;

				int positionOfCurrentExonStartPlusStrand = -1;
				int positionOfCurrentExonStartMinusStrand = -1;

				HashMap<Integer, Integer> mapPositionToCoverageOnPlusStrand = new HashMap<Integer, Integer>();
				HashMap<Integer, Integer> mapPositionToCoverageOnMinusStrand = new HashMap<Integer, Integer>();

				//Loop once per window.
				for (int positionOfWindowEnd = analysisRegion.getCoordinateOfStart(); positionOfWindowEnd <= analysisRegion.getCoordinateOfEnd(); positionOfWindowEnd++) {
					countOfAllAlignmentsProcessedWithMinScore += readInCoverageUpToThisPosition(ermPolyIterator, indexOfCurrentReferenceSequence + 1, positionOfWindowEnd, minAlignmentScore, mapPositionToCoverageOnPlusStrand, mapPositionToCoverageOnMinusStrand, readCounter);

					if (mapPositionToCoverageOnPlusStrand.containsKey(positionOfWindowEnd))
						coverageSumInWindowOnPlusStrand += mapPositionToCoverageOnPlusStrand.get(positionOfWindowEnd);
					if (mapPositionToCoverageOnMinusStrand.containsKey(positionOfWindowEnd))
						coverageSumInWindowOnMinusStrand += mapPositionToCoverageOnMinusStrand.get(positionOfWindowEnd);
					if (positionOfWindowEnd >= widthOfMovingWindow) {

						int positionOldWindowStart = positionOfWindowEnd - widthOfMovingWindow;
						if (mapPositionToCoverageOnPlusStrand.containsKey(positionOldWindowStart))
							coverageSumInWindowOnPlusStrand -= mapPositionToCoverageOnPlusStrand.get(positionOldWindowStart);
						if (mapPositionToCoverageOnMinusStrand.containsKey(positionOldWindowStart))
							coverageSumInWindowOnMinusStrand -= mapPositionToCoverageOnMinusStrand.get(positionOldWindowStart);

						// call exons on plus strand
						if (coverageSumInWindowOnPlusStrand >= (minCoverageToCallAsPartOfExon * widthOfMovingWindow) && positionOfWindowEnd < analysisRegion.getCoordinateOfEnd()) {
							
							if (positionOfCurrentExonStartPlusStrand < 0)
								//Begining of new exon.
								positionOfCurrentExonStartPlusStrand = positionOfWindowEnd - widthOfMovingWindow;
								
						} else if (positionOfCurrentExonStartPlusStrand > -1) {
							//Reached the of an exon region.
							countOfAllPredictedExons++;
							ContiguousGenomicFeature exonPredicted = new ContiguousGenomicFeature(idOfCurrentReferenceSequence, Strand.POSITIVE, positionOfCurrentExonStartPlusStrand, positionOfWindowEnd -1);
							exonPredicted = trimExonEndsByCoverage(exonPredicted, mapPositionToCoverageOnPlusStrand, minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon);
							ContiguousGenomicFeatureList.OverlappingGenomicFeaturesSet exonsThatOverlapThisPrediction = listOfContiguousGenomicFeatures.getSetOfFeaturesOverlapping(exonPredicted, true, minOverlapForCountingPredictedExonAsTruePositive);
							
							double maxOverlapWithAnAnnotatedExonEitherPerspective = Math.max(0.0, Math.max(exonsThatOverlapThisPrediction.getMaxOverlapFractionMyPerspective(), exonsThatOverlapThisPrediction.getMaxOverlapFractionTheirPerspective()));
							boolean predictedExonIsTruePositive = maxOverlapWithAnAnnotatedExonEitherPerspective >= minOverlapForCountingPredictedExonAsTruePositive;
							if (predictedExonIsTruePositive) {
								countOfAllPredictedExonsOverlappingKnownExonOnPlusStrand++;
								sumOfFractionsOfPredictedExonOverlappingAnnotatedExon += exonsThatOverlapThisPrediction.getMaxOverlapFractionTheirPerspective();
							}

							//Iterator<ContiguousGenomicFeature> iteratorOverOverlappingExons = exonsThatOverlapThisPrediction.getFeaturesThatOverlapThisOneWithOverlapAnnotated().iterator();
							StringBuffer exonMatchString = new StringBuffer();
							StringBuffer exonOverlapString = new StringBuffer();
							StringBuffer trOverlapString = new StringBuffer();
							List<ContiguousGenomicFeature> exons = new ArrayList<ContiguousGenomicFeature>(exonsThatOverlapThisPrediction.getFeaturesThatOverlapThisOneWithOverlapAnnotated());
							Collections.sort(exons, ContiguousGenomicFeatureComparator.INSTANCE);
							for (Iterator<ContiguousGenomicFeature> iteratorOverOverlappingExons = exons.iterator(); iteratorOverOverlappingExons.hasNext();) {
								
								ContiguousGenomicFeature exonAnnotated = iteratorOverOverlappingExons.next();
								String keyToExon = idOfCurrentReferenceSequence + "_+_" + exonAnnotated.getCoordinateOfStart() + "_" + exonAnnotated.getCoordinateOfEnd();
								if (mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.containsKey(keyToExon))
									mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.put(keyToExon, Math.max(mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.get(keyToExon), exonAnnotated.getOverlapFractionFromMyPerspective()));
								else
									mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.put(keyToExon, exonAnnotated.getOverlapFractionFromMyPerspective());

								if (exonMatchString.length() > 0) exonMatchString.append(",");
								if (exonOverlapString.length() > 0) exonOverlapString.append(",");
								if (trOverlapString.length() > 0) trOverlapString.append(",");
								exonMatchString.append(exonAnnotated.getCoordinateOfStart()+"-"+exonAnnotated.getCoordinateOfEnd());
								exonOverlapString.append(OVERLAP_FORMAT.format(exonAnnotated.getOverlapFractionFromMyPerspective()));
								trOverlapString.append(OVERLAP_FORMAT.format(exonAnnotated.getOverlapFractionFromItsPerspective()));
							}
							double coverage = calculateCoverage(exonPredicted, mapPositionToCoverageOnPlusStrand);
							writerOutputPlusStrandGFFFile.write(idOfCurrentReferenceSequence + "\t" + "ntr_finder\ttranscript_region\t"
									+ (exonPredicted.getCoordinateOfStart() +1) + "\t" + (exonPredicted.getCoordinateOfEnd() +1) + "\t.\t+\t.\t");
							writerOutputPlusStrandGFFFile.write("ID=PTR_"+countOfAllPredictedExons+"; ");
							writerOutputPlusStrandGFFFile.write("coverage="+COVERAGE_FORMAT.format(coverage)+"; ");
							if (exonMatchString.length() > 0) {
								writerOutputPlusStrandGFFFile.write("atr_match="+exonMatchString+"; ");
								writerOutputPlusStrandGFFFile.write("atr_overlap="+exonOverlapString+"; ");
								writerOutputPlusStrandGFFFile.write("ptr_overlap="+trOverlapString+"; ");
							}
							writerOutputPlusStrandGFFFile.newLine();

							for (int positionCurrent = positionOfCurrentExonStartPlusStrand; positionCurrent <= (positionOfWindowEnd - widthOfMovingWindow); positionCurrent++)
								mapPositionToCoverageOnPlusStrand.remove(positionCurrent);

							positionOfCurrentExonStartPlusStrand = -1;
						} else
							mapPositionToCoverageOnPlusStrand.remove(positionOfWindowEnd - widthOfMovingWindow);


						// call exons on minus strand
						if (coverageSumInWindowOnMinusStrand >= (minCoverageToCallAsPartOfExon * widthOfMovingWindow) && positionOfWindowEnd < analysisRegion.getCoordinateOfEnd()) {

							if (positionOfCurrentExonStartMinusStrand < 0)
								positionOfCurrentExonStartMinusStrand = positionOfWindowEnd - widthOfMovingWindow;
							
						} else if (positionOfCurrentExonStartMinusStrand > -1) {

							countOfAllPredictedExons++;

							ContiguousGenomicFeature exonPredicted = new ContiguousGenomicFeature(idOfCurrentReferenceSequence, Strand.NEGATIVE, positionOfCurrentExonStartMinusStrand, positionOfWindowEnd -1);
							exonPredicted = trimExonEndsByCoverage(exonPredicted, mapPositionToCoverageOnMinusStrand, minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon);
							ContiguousGenomicFeatureList.OverlappingGenomicFeaturesSet exonsThatOverlapThisPrediction = listOfContiguousGenomicFeatures.getSetOfFeaturesOverlapping(exonPredicted, true, minOverlapForCountingPredictedExonAsTruePositive);

							double maxOverlapWithAnAnnotatedExonEitherPerspective = Math.max(0.0, Math.max(exonsThatOverlapThisPrediction.getMaxOverlapFractionMyPerspective(), exonsThatOverlapThisPrediction.getMaxOverlapFractionTheirPerspective()));
							boolean predictedExonIsTruePositive = (maxOverlapWithAnAnnotatedExonEitherPerspective  >= minOverlapForCountingPredictedExonAsTruePositive);
							if (predictedExonIsTruePositive) {
								countOfAllPredictedExonsOverlappingKnownExonOnMinusStrand++;
								sumOfFractionsOfPredictedExonOverlappingAnnotatedExon += exonsThatOverlapThisPrediction.getMaxOverlapFractionTheirPerspective();
							}

							//Iterator<ContiguousGenomicFeature> iteratorOverOverlappingExons = exonsThatOverlapThisPrediction.getFeaturesThatOverlapThisOneWithOverlapAnnotated().iterator();
							StringBuffer exonMatchString = new StringBuffer();
							StringBuffer exonOverlapString = new StringBuffer();
							StringBuffer trOverlapString = new StringBuffer();
							List<ContiguousGenomicFeature> exons = new ArrayList<ContiguousGenomicFeature>(exonsThatOverlapThisPrediction.getFeaturesThatOverlapThisOneWithOverlapAnnotated());
							Collections.sort(exons, ContiguousGenomicFeatureComparator.INSTANCE);
							for (Iterator<ContiguousGenomicFeature> iteratorOverOverlappingExons = exons.iterator(); iteratorOverOverlappingExons.hasNext();) {
								ContiguousGenomicFeature exonAnnotated = iteratorOverOverlappingExons.next();
								String keyToExon = idOfCurrentReferenceSequence + "_-_" + exonAnnotated.getCoordinateOfStart() + "_" + exonAnnotated.getCoordinateOfEnd();
								if (mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.containsKey(keyToExon))
									mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.put(keyToExon, Math.max(mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.get(keyToExon), exonAnnotated.getOverlapFractionFromMyPerspective()));
								else
									mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.put(keyToExon, exonAnnotated.getOverlapFractionFromMyPerspective());
								
								if (exonMatchString.length() > 0) exonMatchString.append(",");
								if (exonOverlapString.length() > 0) exonOverlapString.append(",");
								if (trOverlapString.length() > 0) trOverlapString.append(",");
								exonMatchString.append(exonAnnotated.getCoordinateOfStart()+"-"+exonAnnotated.getCoordinateOfEnd());
								exonOverlapString.append(OVERLAP_FORMAT.format(exonAnnotated.getOverlapFractionFromMyPerspective()));
								trOverlapString.append(OVERLAP_FORMAT.format(exonAnnotated.getOverlapFractionFromItsPerspective()));
							}

							double coverage = calculateCoverage(exonPredicted, mapPositionToCoverageOnMinusStrand);
							writerOutputMinusStrandGFFFile.write(idOfCurrentReferenceSequence + "\t" + "ntr_finder\ttranscript_region\t"
									+ (exonPredicted.getCoordinateOfStart() +1) + "\t" + (exonPredicted.getCoordinateOfEnd() +1) + "\t.\t-\t.\t");
							writerOutputMinusStrandGFFFile.write("ID=PTR_"+countOfAllPredictedExons+"; ");
							writerOutputMinusStrandGFFFile.write("coverage="+COVERAGE_FORMAT.format(coverage)+"; ");
							if (exonMatchString.length() > 0) {
								writerOutputMinusStrandGFFFile.write("atr_match="+exonMatchString+"; ");
								writerOutputMinusStrandGFFFile.write("atr_overlap="+exonOverlapString+"; ");
								writerOutputMinusStrandGFFFile.write("ptr_overlap="+trOverlapString+"; ");
							}
							//+ countOfAllPredictedExons + "_" + new DecimalFormat("0.00").format(maxOverlapWithAnAnnotatedExonEitherPerspective) + "_" + COVERAGE_FORMAT.format(coverage) );
							writerOutputMinusStrandGFFFile.newLine();

							for (int positionCurrent = positionOfCurrentExonStartMinusStrand; positionCurrent <= (positionOfWindowEnd - widthOfMovingWindow); positionCurrent++)
								mapPositionToCoverageOnMinusStrand.remove(positionCurrent);

							positionOfCurrentExonStartMinusStrand = -1;

						} else
							mapPositionToCoverageOnMinusStrand.remove(positionOfWindowEnd - widthOfMovingWindow);

					}
				}
				System.err.println("Number of alignments Processed:\t" + readCounter.getCount());
				System.err.println("Number of alignments processed with min score:\t" + countOfAllAlignmentsProcessedWithMinScore);

			} // End per region loop.

			System.out.println("Number of alignments processed:\t" + readCounter.getCount());
			System.out.println("Number of alignments processed with min score:\t" + countOfAllAlignmentsProcessedWithMinScore);

			double sumOfFractionsOfAnnotatedExonOverlappingPredictedExon = 0.0;
			Iterator<Double> iteratorOverFractionsOverlapped = mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.values().iterator();
			while (iteratorOverFractionsOverlapped.hasNext())
				sumOfFractionsOfAnnotatedExonOverlappingPredictedExon += iteratorOverFractionsOverlapped.next();
			
			double averageFractionOfAnnotatedExonOverlappingPredictedExon= 0;
			if (!mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.isEmpty()) averageFractionOfAnnotatedExonOverlappingPredictedExon = sumOfFractionsOfAnnotatedExonOverlappingPredictedExon / mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.size();

			int countOfAllPredictedExonsOverlappingAnAnnotatedExon = countOfAllPredictedExonsOverlappingKnownExonOnPlusStrand + countOfAllPredictedExonsOverlappingKnownExonOnMinusStrand;
			double averageFractionOfPredictedExonOverlappingAnnotatedExon = 0;
			if (countOfAllPredictedExonsOverlappingAnAnnotatedExon > 0) averageFractionOfPredictedExonOverlappingAnnotatedExon = sumOfFractionsOfPredictedExonOverlappingAnnotatedExon / countOfAllPredictedExonsOverlappingAnAnnotatedExon;

			System.out.println(
					Utilities.join("\t",
					ExonFinder.WIDTH_OF_MOVING_WINDOW,
					ExonFinder.MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON,
					ExonFinder.MIN_OVERLAP_FOR_COUNTING_PREDICTED_EXON_AS_TRUE_POSITIVE,
					ExonFinder.MIN_ALIGNMENT_SCORE,
					ExonFinder.MIN_FRACTION_OF_AVERAGE_COVERAGE_BEFORE_TRIMMING_POSITION_FROM_PREDICTED_EXON,
					ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS,
					ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS,
					ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS_OVERLAPPING_AN_ANNOTATED_EXON,
					ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS,
					ExonFinder.AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON,
					ExonFinder.AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON
					));
			System.out.println(
					Utilities.join("\t",
					widthOfMovingWindow,
					minCoverageToCallAsPartOfExon,
					minOverlapForCountingPredictedExonAsTruePositive,
					minAlignmentScore,
					minFractionOfAverageCoverageBeforeTrimmingPositionFromPredictedExon,
					countOfAllPredictedExons,
					listOfContiguousGenomicFeatures.getSetOfAllFeatures().size(),
					countOfAllPredictedExonsOverlappingAnAnnotatedExon,
					mapAnnotatedExonToMaxOverlapFractionAnnotationPerspective.size(),
					averageFractionOfAnnotatedExonOverlappingPredictedExon,
					averageFractionOfPredictedExonOverlappingAnnotatedExon
					));

		} finally {
			writerOutputPlusStrandGFFFile.close();
			writerOutputMinusStrandGFFFile.close();
		}
	}

	private static int readInCoverageUpToThisPosition(SeekablePolyIterator<ExtendedReadMapping> ermIterator, int indexOfCurrentReferenceSequence, int positionToReadTo, int minAlignmentScore,
			HashMap<Integer, Integer> mapPositionToCoverageOnPlusStrand, HashMap<Integer, Integer> mapPositionToCoverageOnMinusStrand, Counter readCounter) throws IOException {
		int countOfReadsProcessedWithMinAlignScore = 0;
		while (ermIterator.peek() != null &&
				ermIterator.peek().getIndexOfMatchingReferenceSequence() == indexOfCurrentReferenceSequence &&
				Math.abs(ermIterator.peek().getPositionOfAlignmentStartInReferenceSequence()) <= positionToReadTo) {

			ExtendedReadMapping extendedReadMappingCurrent = ermIterator.next();
			readCounter.increment();
			if (showProgress && readCounter.getCount() % 10000 == 0) {
				System.err.printf("%9d Reads Processed", readCounter.getCount());
				System.err.print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			}
			//System.out.println(readCounter.getCount() + " " + extendedReadMappingCurrent);
			if (extendedReadMappingCurrent.getScore() >= minAlignmentScore) {
				countOfReadsProcessedWithMinAlignScore++;
				int positionStartInReferenceAbsolute = Math.abs(extendedReadMappingCurrent.getPositionOfAlignmentStartInReferenceSequence());
				if (extendedReadMappingCurrent.getPositionOfAlignmentStartInReferenceSequence() >= 0) {
					for (int positionInReference = positionStartInReferenceAbsolute; positionInReference < positionStartInReferenceAbsolute + extendedReadMappingCurrent.getLengthOfAlignment(); positionInReference++) {
						if (mapPositionToCoverageOnPlusStrand.containsKey(positionInReference))
							mapPositionToCoverageOnPlusStrand.put(positionInReference, mapPositionToCoverageOnPlusStrand.get(positionInReference) +1);
						else
							mapPositionToCoverageOnPlusStrand.put(positionInReference, 1);
					}
				} else {
					for (int positionInReference = positionStartInReferenceAbsolute; positionInReference < positionStartInReferenceAbsolute + extendedReadMappingCurrent.getLengthOfAlignment(); positionInReference++) {
						if (mapPositionToCoverageOnMinusStrand.containsKey(positionInReference))
							mapPositionToCoverageOnMinusStrand.put(positionInReference, mapPositionToCoverageOnMinusStrand.get(positionInReference) +1);
						else
							mapPositionToCoverageOnMinusStrand.put(positionInReference, 1);
					}
				}
			}
		}
		return countOfReadsProcessedWithMinAlignScore;
	}

	private static ContiguousGenomicFeature trimExonEndsByCoverage(ContiguousGenomicFeature exon, HashMap<Integer, Integer> mapPositionToCoverage, double minFractionOfAverageCoverageBeforeTrimmingPositionFromExon) throws Exception {
		ContiguousGenomicFeature exonTrimmed = (ContiguousGenomicFeature)exon.clone();    //new ContiguousGenomicFeature(exon.getCoordinateOfStart(), exon.getCoordinateOfEnd());
		double coverageInExon = 0.0;
		for (int positionCurrent = exon.getCoordinateOfStart(); positionCurrent <= exon.getCoordinateOfEnd(); positionCurrent++)
			if (mapPositionToCoverage.containsKey(positionCurrent))
				coverageInExon += mapPositionToCoverage.get(positionCurrent);
		double averageCoverageAcrossExon = coverageInExon / exon.getLength();

		int positionCurrent = exon.getCoordinateOfStart();
		int coverageAtCurrentPosition = 0;
		if (mapPositionToCoverage.containsKey(positionCurrent))
			coverageAtCurrentPosition = mapPositionToCoverage.get(positionCurrent);
		while (positionCurrent <= exon.getCoordinateOfEnd()
				&& (coverageAtCurrentPosition / averageCoverageAcrossExon) < minFractionOfAverageCoverageBeforeTrimmingPositionFromExon) {
			positionCurrent++;
			if (mapPositionToCoverage.containsKey(positionCurrent))
				coverageAtCurrentPosition = mapPositionToCoverage.get(positionCurrent);
			else
				coverageAtCurrentPosition = 0;
		}
		exonTrimmed.setCoordinateOfStart(positionCurrent);
		//System.out.println(exonTrimmed.getCoordinateOfStart() - exon.getCoordinateOfStart() + " trimmed from 5' end");

		positionCurrent = exon.getCoordinateOfEnd();
		coverageAtCurrentPosition = 0;
		if (mapPositionToCoverage.containsKey(positionCurrent))
			coverageAtCurrentPosition = mapPositionToCoverage.get(positionCurrent);
		while (positionCurrent >= exonTrimmed.getCoordinateOfStart()
				&& (coverageAtCurrentPosition / averageCoverageAcrossExon) < minFractionOfAverageCoverageBeforeTrimmingPositionFromExon) {
			positionCurrent--;
			if (mapPositionToCoverage.containsKey(positionCurrent))
				coverageAtCurrentPosition = mapPositionToCoverage.get(positionCurrent);
			else
				coverageAtCurrentPosition = 0;
		}
		exonTrimmed.setCoordinateOfEnd(positionCurrent);
		//System.out.println(exon.getCoordinateOfEnd() - exonTrimmed.getCoordinateOfEnd() + " trimmed from 3' end");

		return exonTrimmed;
	}


	public static HashMap<Integer, int[]> loadFeatureStartToEndsFromGFF(File fileFeatureAnnotationGFF, String nameOfChromosome, char strand, boolean mapEndToStart) throws IOException {

		HashMap<Integer, int[]> mapCoordinateToCoordinates = new HashMap<Integer, int[]>();  // either startToEnd or EndToStart depending on value of mapEndToStart

		BufferedReader readerFaetureAnnotationFile = new BufferedReader(new FileReader(fileFeatureAnnotationGFF));
		String headerOfRead;
		while ((headerOfRead = readerFaetureAnnotationFile.readLine()) != null) {
			if (!headerOfRead.startsWith("#")) {
				String tokens[] = headerOfRead.split("\t");
				if (tokens[0].equalsIgnoreCase(nameOfChromosome) && tokens[6].equalsIgnoreCase(strand + "")) {

					int coordinateStart = Integer.parseInt(tokens[3]);
					int coordinateEnd = Integer.parseInt(tokens[4]);

					if (mapEndToStart) {

						if (mapCoordinateToCoordinates.containsKey(coordinateEnd)) {
							int starts[] = mapCoordinateToCoordinates.get(coordinateEnd);
							int startsNew[] = Arrays.copyOf(starts, starts.length +1);
							startsNew[starts.length] = coordinateStart;
							mapCoordinateToCoordinates.put(coordinateEnd, startsNew);
						} else {
							int startsNew[] = new int[1];
							startsNew[0] = coordinateStart;
							mapCoordinateToCoordinates.put(coordinateEnd, startsNew);
						}

					} else {

						if (mapCoordinateToCoordinates.containsKey(coordinateStart)) {
							int ends[] = mapCoordinateToCoordinates.get(coordinateStart);
							int endsNew[] = Arrays.copyOf(ends, ends.length +1);
							endsNew[ends.length] = coordinateEnd;
							mapCoordinateToCoordinates.put(coordinateStart, endsNew);
						} else {
							int endsNew[] = new int[1];
							endsNew[0] = coordinateEnd;
							mapCoordinateToCoordinates.put(coordinateStart, endsNew);
						}

					}
				}
			}
		}
		readerFaetureAnnotationFile.close();

		return mapCoordinateToCoordinates;

	}
	
	private static float calculateCoverage(ContiguousGenomicFeature feature, Map<Integer, Integer> coverageMap) {
		int start = feature.getCoordinateOfStart();
		int end = feature.getCoordinateOfEnd();
		if (start > end) {
			int tmp = start;
			start = end;
			end = tmp;
		}
		int baseCount = 0;
		for (int i=start; i<=end; i++) {
			Integer cov = coverageMap.get(i);
			if (cov != null) baseCount += cov;
		}
		return (float)baseCount / feature.getLength();
	}
	
}

/** substitute for an int pointer */
class Counter {
	private int count = 0;
	public void increment() { count++; }
	public int getCount() { return count; } 
}

class FeatureKey {
	
	private ContiguousGenomicFeature feature;
	
	public FeatureKey(ContiguousGenomicFeature feature ) {
		this.feature = feature;
	}
	
	@Override
	public String toString() {
		return feature.getIdOfReferenceSequence()+feature.getCoordinateOfStart()+":"+feature.getCoordinateOfEnd()+feature.getStrand();
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		return toString().equals(obj);
	}
}
