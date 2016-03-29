package com.lifetechnologies.solid.wt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.config.ConfigKey;

/**
 * Launches NTR finding jobs and merges the results.
 * @author mullermw
 *
 */
public class NTRFinderLauncher {

	private static Logger logger = Logger.getLogger(NTRFinderLauncher.class.toString());
	
	public void launch(Map<ConfigKey, Set<String>> config) throws Exception {
		this.launch(config, true);
	}
	
	/**
	 * Launch job using the specified config.  config should be validated before
	 * calling launch().
	 * @param config
	 * @param submitJobs if true, jobs will be sent to the cluster.
	 * @throws Exception
	 */
	public void launch(Map<ConfigKey, Set<String>> config, boolean submitJobs) throws Exception {
		File fileReferenceFasta = new File(Utilities.firstValue(config.get(ConfigKey.WT_FILE_REFERENCE)));
		File outputFolder = new File(Utilities.firstValue(config.get(ConfigKey.WT_OUTPUT_DIR)));
		if (!outputFolder.exists()) outputFolder.mkdir();
		File tmpFolder = new File(outputFolder, "tmp");
		if (!tmpFolder.exists()) tmpFolder.mkdir();
		
		List<ContiguousGenomicFeature> regionsToAnalyze = new ArrayList<ContiguousGenomicFeature>();
		if (config.containsKey(ConfigKey.WT_NTR_GENOMIC_REGION))
			for (String str : config.get(ConfigKey.WT_NTR_GENOMIC_REGION))
				regionsToAnalyze.addAll(ContiguousGenomicFeature.parseFeatures(str));
				
		File referenceIndexFile = new File(tmpFolder, fileReferenceFasta.getName().concat(".idx"));
		IndexedFastaDatabase fastaDatabaseOfReference = new IndexedFastaDatabase(fileReferenceFasta, referenceIndexFile);
		SortedMap<String, Long> mapHeaderToSequenceLength = fastaDatabaseOfReference.getSequenceLengths();
		Map<String, Long> idToSequenceLengths = new HashMap<String, Long>();
		for (Map.Entry<String, Long> entry: mapHeaderToSequenceLength.entrySet())
			idToSequenceLengths.put(entry.getKey().replaceAll("\\s+.*", ""), entry.getValue());
		
		for (ContiguousGenomicFeature region : regionsToAnalyze) {
			if (region.getCoordinateOfStart() < 0 ) region.setCoordinateOfStart(0);
			if (region.getCoordinateOfEnd() < 0) region.setCoordinateOfEnd(idToSequenceLengths.get(region.getIdOfReferenceSequence()).intValue() - 1);
		}
		
		if (regionsToAnalyze.isEmpty()) {
			List<String> listOfHeaders = fastaDatabaseOfReference.getListOfHeaders(true, false);
			Set<String> ids = new HashSet<String>();
			for (String header : listOfHeaders) {
				String id = header.replaceAll("\\s+.*", "");
				if (ids.contains(id)) {
					logger.warning("Ignoring duplicate sequence id: '"+id+"'.");
				} else {
					regionsToAnalyze.add(new ContiguousGenomicFeature(id, Strand.EITHER, 0, mapHeaderToSequenceLength.get(header).intValue() - 1));
					ids.add(id);
				} 
			}
		}
		Collections.sort(regionsToAnalyze, ContiguousGenomicFeatureComparator.INSTANCE);
		
		long sizeOfReference = 0;
		for (ContiguousGenomicFeature region : regionsToAnalyze)
			sizeOfReference += region.getLength();
		Integer maxNTRsPerMBase = DecimalFormat.getIntegerInstance().parse(Utilities.firstValue(config.get(ConfigKey.WT_NTR_MAX_PTRS_PER_MEGABASE))).intValue();
		Integer maxNTRs = Math.round(maxNTRsPerMBase * (sizeOfReference / 1E6f));

		//Write the names of the sorted max files to a file.
        File sortedMaxListFile = new File(tmpFolder, "sorted_max_files");
        PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(sortedMaxListFile)));
        try {
        	for (String sortedMaxFile : config.get(ConfigKey.WT_NTR_MAX_FILE))
        		writer.println(new File(sortedMaxFile).getAbsolutePath());
        } finally {
        	writer.close();
        }
		
        //Setup the cluster submission
		JobSubmissionParameters params = new JobSubmissionParameters();
        params.setEnvironment(Utilities.firstValue(config.get(ConfigKey.QUEUE_SYS)));
        params.setQueueName(Utilities.firstValue(config.get(ConfigKey.QUEUE_SYS_QUEUE)));
        params.setResourceString(Utilities.firstValue(config.get(ConfigKey.QUEUE_SYS_RESOURCE_STRING)));
        params.setAdditionalOptions(Utilities.firstValue(config.get(ConfigKey.QUEUE_SYS_OPTIONS)));
        params.setMemoryRequirement(new Double(1.5e9).longValue());
        params.setRerunnable(false);
        ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
        
        File wtHome = new File (System.getProperty("com.lifetechnologies.solid.wt.home"));
        File libDir = new File (wtHome, "lib");
        File pkgDir = new File (wtHome, "pkg");
		File refFile = new File(Utilities.firstValue(config.get(ConfigKey.WT_FILE_REFERENCE)));
		File exonRefFile = null;
		if (config.get(ConfigKey.WT_NTR_FILE_ATR_REFERENCE) != null)
			exonRefFile = new File(Utilities.firstValue(config.get(ConfigKey.WT_NTR_FILE_ATR_REFERENCE)));
		
        String classpath = "'"+libDir.getAbsolutePath()+"/*:"+pkgDir.getAbsolutePath()+"/*'";
        String memoryRequirement = (int)Math.ceil((double)clusterInterface.getJobSubmissionParameters().getMemoryRequirement() / Constants.BYTES_PER_MEGABYTE) + "m";
        
		//Each config in runConfigs is a separate job.
		List<File> jobOutputFiles = new ArrayList<File>();

		Double overlap = Double.valueOf(Utilities.firstValue(config.get(ConfigKey.WT_NTR_MIN_OVERLAP)));
		Integer minScore = Integer.valueOf(Utilities.firstValue(config.get(ConfigKey.WT_NTR_MIN_ALIGNMENT_SCORE)));
		Double trimmingFraction = Double.valueOf(Utilities.firstValue(config.get(ConfigKey.WT_NTR_TRIMMING_FRACTION)));
		List<Integer> windowSizes = Utilities.toIntegerList(config.get(ConfigKey.WT_NTR_MIN_WINDOW_SIZE));
		Collections.sort(windowSizes);
		List<Double> minWindowCoverages = Utilities.toDoubleList(config.get(ConfigKey.WT_NTR_MIN_WINDOW_COVERAGE));
		Collections.sort(minWindowCoverages);
		
		logger.info("Sumbitting Jobs");
		for (ContiguousGenomicFeature regionToAnalyze : regionsToAnalyze) {
			for (Integer windowSize : windowSizes) {
				for (Double minWindowCoverage : minWindowCoverages) {
					String baseName = Utilities.join("_", regionToAnalyze.toBriefString(), windowSize, minWindowCoverage, overlap, minScore, trimmingFraction);
					baseName = Utilities.toStringSafeForFile(baseName);
					File resultDir = new File(tmpFolder, baseName);
					resultDir.mkdirs();
					File summaryFile = new File(resultDir, baseName+".summary");
					File scriptFile = new File(resultDir, baseName+".sh").getAbsoluteFile();
					String cmdStr = String.format("java -cp %s -Xmx%s %s %s %s %s %s %s %s %s %s %s %s %s %s > %s",
			        	classpath,
			        	memoryRequirement,
			        	//We want all jobs using the same index file.
			        	"-D"+ExonFinder.INDEX_DIR_SYS_PROPERTY+"="+referenceIndexFile.getParentFile().getAbsolutePath(),
			        	"-D"+ExonFinder.SHOW_PROGRESS_SYS_PROPERTY+"=false",
			        	ExonFinder.class.getName(), 
			        	windowSize,
			        	minWindowCoverage,
			        	overlap,
			        	minScore,
			        	trimmingFraction,
			        	sortedMaxListFile.getAbsolutePath(),
			        	refFile.getAbsolutePath(),
			        	exonRefFile == null ? "\"\"" : exonRefFile.getAbsolutePath(),
			        	resultDir.getAbsolutePath(),
			        	regionToAnalyze.toBriefString(),
			        	summaryFile.getAbsolutePath());
					ArrayList<String> cmd = new ArrayList<String>();
					cmd.add(cmdStr);
					logger.fine("Submitting job for NTR_WINDOW_SIZE="+windowSize+" NTR_MIN_WINDOW_COVERAGE"+minWindowCoverage+" NTR_GENOMIC_REGION="+regionToAnalyze.toBriefString()+"\n"+scriptFile);
					if (submitJobs) jobOutputFiles.add(clusterInterface.executeJob(scriptFile, cmd));
				}
			}
        }
		
		if (submitJobs) {
			File fileJobRemoval = new File(outputFolder + "/remove_jobs.sh");
			clusterInterface.writeMasterJobRemovalFileFromLog(fileJobRemoval);
			fileJobRemoval.setExecutable(true);
	        logger.info("Waiting for jobs to complete");
	        while (!clusterInterface.checkIfLoggedJobsComplete()) {
	        	Thread.sleep(5000);
	        }
        }
        
        List<File> failures = clusterInterface.getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure();
        if (failures.isEmpty()) {
        	logger.info("All jobs completed successfully.");
        } else {
        	StringBuffer msg = new StringBuffer("Some of your jobs failed.  Check these files:\n");
        	for (File file : failures)
        		msg.append(                     "  "+file+"\n");
        	if (failures.size() < jobOutputFiles.size())
        		msg.append(                     "The NTR finding process is proceeding with succesful results.");
        	logger.warning(msg.toString());
        	if (failures.size() >= jobOutputFiles.size()) {
        		throw new Exception("All Jobs failed.  Cannot continue.");
        	}
        }

        logger.info("Merging NTR reports");
        //Merge reports into final result.
		List<Map<String, Number>> summaries = new ArrayList<Map<String, Number>>();
		int numFilesToMerge = windowSizes.size() * minWindowCoverages.size() * regionsToAnalyze.size();
		int filesMerged = 0;
    	for (Integer windowSize : windowSizes) {
    		for (Double minWindowCoverage : minWindowCoverages) {
        		Map<String, Number> mergedSummary = null;
        		for (ContiguousGenomicFeature regionToAnalyze : regionsToAnalyze) {
        			filesMerged++;
        			logger.fine("Summary Merge #"+filesMerged+" of "+numFilesToMerge+"): window size=" + windowSize + ", min window coverage="+minWindowCoverage);
        			String baseName = Utilities.join("_", regionToAnalyze.toBriefString(), windowSize, minWindowCoverage, overlap, minScore, trimmingFraction);
        			baseName = Utilities.toStringSafeForFile(baseName);
        			File resultDir = new File(tmpFolder, baseName);
					File summaryFile = new File(resultDir, baseName+".summary");
					
					logger.fine("Parsing summary file "+summaryFile);
					try {
						HashMap<String, Number> thisSummary = parseSummaryFile(summaryFile);
						if (mergedSummary == null) {
							mergedSummary = clone(thisSummary);
	        			} else {
	        				mergeSummary(thisSummary, mergedSummary);
	        			}
					} catch (Exception e) {
						logger.severe("Failed to processs summary file: "+summaryFile);
					}
        		}
        		summaries.add(mergedSummary);
        	}
        }
    	
    	logger.info("Writing report.");
        List<String> reportKeys = Arrays.asList(new String[] {
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
        });
        File masterSummaryFile = new File(outputFolder, "NTR_report.txt");
        PrintWriter out = new PrintWriter(new FileOutputStream(masterSummaryFile));
        try {
        	out.println(Utilities.join(reportKeys, "\t"));
        	for (Map<String, Number> summary : summaries) {
        		for (Iterator<String> it=reportKeys.iterator(); it.hasNext();) {
        			String key = it.next();
        			Number value = summary.get(key);
        			out.print(value);
        			if (it.hasNext()) out.print("\t");
        		}
        		out.println();
        	}
        } finally {
        	out.close();
        }
    	
        logger.info("Merging Gff results");
        //Merge job results into final result.
		filesMerged = 0;
		File resultsDir = new File(outputFolder, "ptr_results");
		for (Map<String, Number> summary : summaries) {
			Integer windowSize = summary.get(ExonFinder.WIDTH_OF_MOVING_WINDOW).intValue();
			Double minWindowCoverage = summary.get(ExonFinder.MIN_COVERAGE_TO_CALL_AS_PART_OF_EXON).doubleValue();
			Integer numNTRs = summary.get(ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS).intValue();
			// Too many NTRs found.  Don't copy the GFF.
			if (numNTRs > maxNTRs) {
				filesMerged+= regionsToAnalyze.size();
				logger.info(numNTRs + " > "+maxNTRs+" NTRs found for window size = " + windowSize + ", min window coverage="+minWindowCoverage+".  Skipping GFF Merge.");
				continue;
			}
			File windowSizeDir = new File(resultsDir, "window_size_"+windowSize);
			File windowCoverageDir = new File(windowSizeDir, "window_coverage_"+minWindowCoverage);
			windowCoverageDir.mkdirs();
			File mergedPlusGffFile = new File(windowCoverageDir, Utilities.join("_", "NTRFinder."+windowSize, minWindowCoverage, overlap, minScore, trimmingFraction+".plus.gff"));
			File mergedMinusGffFile = new File(mergedPlusGffFile.getParent(), mergedPlusGffFile.getName().replace(".plus.gff", ".minus.gff"));		
			PrintWriter plusWriter = new PrintWriter(new BufferedWriter(new FileWriter(mergedPlusGffFile)));
			plusWriter.printf("track name=NTRFinderPlus description=\"NTRFinder (+) window-size_%d window-coverage_%1.1f min-score_%d trimming-fraction_%1.2f\"\n", windowSize, minWindowCoverage, minScore, trimmingFraction);
			PrintWriter minusWriter = new PrintWriter(new BufferedWriter(new FileWriter(mergedMinusGffFile)));
			minusWriter.printf("track name=NTRFinderMinus description=\"NTRFinder (-) window-size_%d window-coverage_%1.1f min-score_%d trimming-fraction_%1.2f\"\n", windowSize, minWindowCoverage, minScore, trimmingFraction);
			try {
				for (ContiguousGenomicFeature regionToAnalyze : regionsToAnalyze) {
					filesMerged++;
					logger.fine("GFF Merge #"+filesMerged+" of "+numFilesToMerge+"): window size=" + windowSize + ", min window coverage="+minWindowCoverage);
					String baseName = Utilities.join("_", regionToAnalyze.toBriefString(), windowSize, minWindowCoverage, overlap, minScore, trimmingFraction);
					baseName = Utilities.toStringSafeForFile(baseName);
					File resultDir = new File(tmpFolder, baseName);

					File plusGffFile = new File(resultDir, Utilities.join("_", "exonFinder."+windowSize, minWindowCoverage, overlap, minScore, trimmingFraction.toString().concat(".plus.gff")));
					logger.fine("Concatenating gff file: "+plusGffFile);
					int ptrCounter = 1;
					BufferedReader plusReader = new BufferedReader(new FileReader(plusGffFile));
					try {
						String line;
						//Skip over headers.
						for (line = plusReader.readLine(); line != null; line = plusReader.readLine())
							if (line.split("\t").length >= 9) break;

						for (; line != null; line = plusReader.readLine()) {
							String[] fields = line.split("\t");
							fields[8] = fields[8].replaceFirst("ID=.*?;", "ID=PTR_"+(ptrCounter++)+";");
							line = Utilities.join("\t", fields);
							plusWriter.println(line);
						}

					} finally {
						plusReader.close();
					}

					File minusGffFile = new File(resultDir, plusGffFile.getName().replace(".plus.gff", ".minus.gff"));
					logger.fine("Concatenating gff file: "+minusGffFile);
					BufferedReader minusReader = new BufferedReader(new FileReader(minusGffFile));
					try {
						String line;
						for (line = minusReader.readLine(); line != null; line = minusReader.readLine())
							if (line.split("\t").length >= 9) break;
						for (; line != null; line = minusReader.readLine()) {
							String[] fields = line.split("\t");
							fields[8] = fields[8].replaceFirst("ID=.*?;", "ID=PTR_"+(ptrCounter++)+";");
							line = Utilities.join("\t", fields);
							minusWriter.println(line);
						}
					} finally {
						minusReader.close();
					}
				}
			} finally {
				plusWriter.close();
				minusWriter.close();
			}
		}
        
        logger.info("Creating plots.");
        createPlots(summaries, outputFolder, maxNTRs);
        
        if (Boolean.valueOf(Utilities.firstValue(config.get(ConfigKey.WT_DELETE_INTERMEDIATE_FILES))))
        	Utilities.deleteDirectory(tmpFolder);
        
	}
		
	/**
	 * return a copy of map.
	 * @param map
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private static Map<String, Number> clone(HashMap<String, Number> map) {
		return (Map<String, Number>)map.clone();
	}
	
	/**
	 * parse an NTR Finding Summary File (Defined in ExonFinder).
	 * @param summaryFile
	 * @return the parsed result as a map.
	 * @throws IOException
	 * @throws ParseException
	 */
	private HashMap<String, Number> parseSummaryFile(File summaryFile) throws IOException, ParseException {
		HashMap<String, Number> result = new HashMap<String, Number>();
		BufferedReader reader = new BufferedReader(new FileReader(summaryFile));
		try {
			List<String> keys = null;
			NumberFormat format = new DecimalFormat();
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				//line with 3 or more tab delimited numbers.
				if (line.matches("^(?:[\\d\\.]+\\t){2,}[\\d\\.]+$")) {
					if (keys != null) {
						String[] values = line.split("\\t");
						Iterator<String> it = keys.iterator();
						for (String valueStr : values) {
							Number value = format.parse(valueStr);
							String key = it.next();
							if (key != null) {
								result.put(key, value);
							}
						}
					}
					keys = null;
				//line with 3 or more tab delimited strings.
				} else if (line.matches("^(?:[^\\t]+\\t){2,}[^\\t]+$")) {
					keys = Arrays.asList(line.split("\\t"));
				}
			}
		} finally {
			reader.close();
		}
		return result;
		
	}
	
	/**
	 * Combine the values of src into dest.  Counts are added up.  Averages are properly merged.
	 * @param src
	 * @param dest
	 */
	private static void mergeSummary(Map<String, Number> src, Map<String, Number> dest) {
			
		//Merging averageFractionOfAnnotatedExonOverlappingPredictedExon
		double currAvgFracAOverP = dest.get(ExonFinder.AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON).doubleValue();
		double srcAvgFracAOverP = src.get(ExonFinder.AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON).doubleValue();
		double currNumAOverP = dest.get(ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS).intValue();
		double srcNumAOverP = src.get(ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS).intValue();
		double currFracSumAOverP = currAvgFracAOverP * currNumAOverP;
		double srcFracSumAOverP = srcAvgFracAOverP * srcNumAOverP;
		double newAvgFracAOverP = 0;
		if (currNumAOverP + srcNumAOverP > 0) newAvgFracAOverP = (currFracSumAOverP + srcFracSumAOverP) / (currNumAOverP + srcNumAOverP);
		dest.put(ExonFinder.AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON, newAvgFracAOverP);
		
		//Merging averageFractionOfPredictedExonOverlappingAnnotatedExon
		double currAvgFracPOverA = dest.get(ExonFinder.AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON).doubleValue();
		double srcAvgFracPOverA = src.get(ExonFinder.AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON).doubleValue();
		double currNumPOverA = dest.get(ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS_OVERLAPPING_AN_ANNOTATED_EXON).intValue();
		double srcNumPOverA = src.get(ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS_OVERLAPPING_AN_ANNOTATED_EXON).intValue();
		double currFracSumPOverA = currAvgFracPOverA * currNumPOverA;
		double srcFracSumPOverA = srcAvgFracPOverA * srcNumPOverA;
		double newAvgFracPOverA = 0;
		if (currNumPOverA + srcNumPOverA > 0) newAvgFracPOverA = (currFracSumPOverA + srcFracSumPOverA) / (currNumPOverA + srcNumPOverA);
		dest.put(ExonFinder.AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON, newAvgFracPOverA);
		
		List<String> addupIntKeys = Arrays.asList(new String[] {
				ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS,
				ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS,        
				ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS_OVERLAPPING_AN_ANNOTATED_EXON,      
				ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS
		});
		for (String key : addupIntKeys) {
			int currValue = dest.get(key).intValue();
			int srcValue = src.get(key).intValue();
			dest.put(key, currValue + srcValue);
		}
	}
	
	/**
	 * Create NTR Finding optimization plots.
	 * @param summaries
	 * @param dest
	 */
	private static void createPlots(Collection<Map<String, Number>> summaries, File dest, int maxNTRsShown) {
		Map<Integer, XYSeries> windowSize2annotFoundSeries = new TreeMap<Integer, XYSeries>();
		Map<Integer, XYSeries> windowSize2POverASeries = new TreeMap<Integer, XYSeries>();
		Map<Integer, XYSeries> windowSize2AOverPSeries = new TreeMap<Integer, XYSeries>();
		double minFractionAnnotatedExonFound = Double.MAX_VALUE;
		double minAverageFractionOfPredictedExonOverlappingAnnotatedExon = Double.MAX_VALUE;
		double minAverageFractionOfAnnotatedExonOverlappingPredictedExon = Double.MAX_VALUE;
		Integer maxNumNTRs = 0;
		for (Map<String, Number> summary : summaries) {
			int windowSize = summary.get(ExonFinder.WIDTH_OF_MOVING_WINDOW).intValue();
			int countOfAllPredictedExons = summary.get(ExonFinder.COUNT_OF_ALL_PREDICTED_EXONS).intValue();
			maxNumNTRs = countOfAllPredictedExons > maxNumNTRs ? countOfAllPredictedExons : maxNumNTRs;
			int countOfAllAnnotatedExons = summary.get(ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS).intValue();
			int countOfAllAnnotatedExonsOverlappingPredictedExons = summary.get(ExonFinder.COUNT_OF_ALL_ANNOTATED_EXONS_OVERLAPPING_PREDICTED_EXONS).intValue();
			double averageFractionOfAnnotatedExonOverlappingPredictedExon = summary.get(ExonFinder.AVERAGE_FRACTION_OF_ANNOTATED_EXON_OVERLAPPING_PREDICTED_EXON).doubleValue();
			double averageFractionOfPredictedExonOverlappingAnnotatedExon = summary.get(ExonFinder.AVERAGE_FRACTION_OF_PREDICTED_EXON_OVERLAPPING_ANNOTATED_EXON).doubleValue();
			
			double fracAnnotedExonFound = 0;
			if (countOfAllAnnotatedExons > 0) fracAnnotedExonFound = (double)countOfAllAnnotatedExonsOverlappingPredictedExons / countOfAllAnnotatedExons;
			
			minFractionAnnotatedExonFound = fracAnnotedExonFound < minFractionAnnotatedExonFound ? fracAnnotedExonFound : minFractionAnnotatedExonFound;
			minAverageFractionOfPredictedExonOverlappingAnnotatedExon = averageFractionOfPredictedExonOverlappingAnnotatedExon < minAverageFractionOfPredictedExonOverlappingAnnotatedExon ? averageFractionOfPredictedExonOverlappingAnnotatedExon : minAverageFractionOfPredictedExonOverlappingAnnotatedExon;
			minAverageFractionOfAnnotatedExonOverlappingPredictedExon = averageFractionOfAnnotatedExonOverlappingPredictedExon < minAverageFractionOfAnnotatedExonOverlappingPredictedExon ? averageFractionOfAnnotatedExonOverlappingPredictedExon : minAverageFractionOfAnnotatedExonOverlappingPredictedExon;
			
			XYSeries annotFoundSeries = windowSize2annotFoundSeries.get(windowSize);
			if (annotFoundSeries == null) {
				if (windowSize2annotFoundSeries.isEmpty())
					annotFoundSeries = new XYSeries("Window Size="+windowSize);
				else
					annotFoundSeries = new XYSeries(windowSize);
				windowSize2annotFoundSeries.put(windowSize, annotFoundSeries);
			}
			annotFoundSeries.add(countOfAllPredictedExons, fracAnnotedExonFound * 100);
			
			XYSeries pOverASeries = windowSize2POverASeries.get(windowSize);
			if (pOverASeries == null) {
				if (windowSize2POverASeries.isEmpty())
					pOverASeries = new XYSeries("Window Size="+windowSize);
				else
					pOverASeries = new XYSeries(windowSize);
				windowSize2POverASeries.put(windowSize, pOverASeries);
			}
			pOverASeries.add(countOfAllPredictedExons, averageFractionOfPredictedExonOverlappingAnnotatedExon);
			
			XYSeries aOverPSeries = windowSize2AOverPSeries.get(windowSize);
			if (aOverPSeries == null) {
				if (windowSize2AOverPSeries.isEmpty())
					aOverPSeries = new XYSeries("Window Size="+windowSize);
				else
					aOverPSeries = new XYSeries(windowSize);
				windowSize2AOverPSeries.put(windowSize, aOverPSeries );
			}
			aOverPSeries.add(countOfAllPredictedExons, averageFractionOfAnnotatedExonOverlappingPredictedExon);
		}
		
		XYSeriesCollection annotFoundData = new XYSeriesCollection();
		for (XYSeries series : windowSize2annotFoundSeries.values())
			annotFoundData.addSeries(series);
		
		XYSeriesCollection pOverAData = new XYSeriesCollection();
		for (XYSeries series : windowSize2POverASeries.values())
			pOverAData.addSeries(series);
		
		XYSeriesCollection aOverPData = new XYSeriesCollection();
		for (XYSeries series : windowSize2AOverPSeries.values())
			aOverPData.addSeries(series);
		
		JFreeChart annotFoundChart = ChartFactory.createXYLineChart("Number of PTRs vs. Fraction of ATRs found", 
				"Number of PTRs",
				"% of ATRs Found",
				annotFoundData, 
				PlotOrientation.VERTICAL, true, false, false);
		
		JFreeChart pOverAChart = ChartFactory.createXYLineChart("Number of PTRs vs. Overlap", 
				"Number of PTRs", 
				"Average fraction of PTR overlapping with ATR", 
				pOverAData, PlotOrientation.VERTICAL, true, false, false);
		
		JFreeChart aOverPChart = ChartFactory.createXYLineChart("Number of TRs vs. Overlap", 
				"Number of PTRs", 
				"Average fraction of ATR overlapping with PTR", 
				aOverPData, PlotOrientation.VERTICAL, true, false, false);

       
        for (JFreeChart chart : new JFreeChart[] {annotFoundChart, pOverAChart, aOverPChart}) {
        	XYPlot plot = chart.getXYPlot();
        	plot.setRenderer(new XYLineAndShapeRenderer(true, true));
        	chart.setBackgroundPaint(plot.getBackgroundPaint());
        	chart.setAntiAlias(true);
        	Axis rangeAxis = plot.getRangeAxis();
        	rangeAxis.setLabelFont(rangeAxis.getLabelFont().deriveFont(10f));
        	if (maxNumNTRs > maxNTRsShown)
        		plot.getDomainAxis().setRange(0, maxNTRsShown);
        }
        
        annotFoundChart.getXYPlot().getRangeAxis().setRange((minFractionAnnotatedExonFound - 0.05)*100, 100);
        pOverAChart.getXYPlot().getRangeAxis().setRange(minAverageFractionOfPredictedExonOverlappingAnnotatedExon - 0.05, 1.0);
        aOverPChart.getXYPlot().getRangeAxis().setRange(minAverageFractionOfAnnotatedExonOverlappingPredictedExon - 0.05, 1.0);
        
        GraphingUtilities.convertToPdf(annotFoundChart, 600, 400, new File(dest, "plot_atrs_found.pdf"));
        GraphingUtilities.convertToPdf(pOverAChart, 600, 400, new File(dest, "plot_ptr_over_atr.pdf"));
        GraphingUtilities.convertToPdf(aOverPChart, 600, 400, new File(dest, "plot_atr_over_ptr.pdf"));
	}
}

