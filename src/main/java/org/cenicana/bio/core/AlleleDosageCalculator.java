package org.cenicana.bio.core;
import org.cenicana.bio.utils.FileUtils;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.concurrent.*;
import java.util.PriorityQueue;

/**
 * Calculates allele dosage matrices from Variant Call Format (VCF) files.
 * 
 * Biological Purpose:
 * In polyploid organisms (e.g., sugarcane which has a ploidy of 10), individuals can 
 * carry multiple copies of a given allele. A simple heterozygous/homozygous 
 * call is insufficient. Instead, we compute the "dosage", which is the proportion 
 * of the chromosomes that carry the reference allele.
 * 
 * Mathematical Approach:
 * Dosage = reads_reference / (reads_reference + reads_alternate)
 * 
 * Phase 1 features a Statistical Depth Filter: 
 * If the total number of reads (Ref + Alt) is below a user-defined threshold, the 
 * genotype is considered statistically unreliable (too much noise) and is treated 
 * as a missing value, to be subsequently imputed based on population data.
 * 
 * Phase 2 features K-Nearest Neighbors (KNN) Imputation:
 * Instead of imputing missing values with the global mode/mean, it finds the K 
 * individuals most genetically similar to the sample with the missing data, 
 * and imputes using their average dosage. This is implemented via a low-memory 
 * two-pass algorithm to support massive VCF files.
 * 
 * Phase 3 features Adaptive Rounding via 1D K-Means Clustering:
 * Instead of strictly rounding raw frequencies to the nearest ploidy level (which 
 * suffers from reference bias), this uses Unsupervised Machine Learning to find 
 * the true density clusters of the data, assigning individuals to their natural group.
 */
public class AlleleDosageCalculator {

	private int numSNPs = 0;
	private int numGenotypes = 0;
	public int maxSnps = Integer.MAX_VALUE;
	public java.util.Set<String> subsetSnps = null;

	public static class DosageResult {
		public String snpId;
		public double[] dosages;
		public int[] refDepths;
		public int[] altDepths;
	}
	
	// Formatter to round decimal representations to 6 decimal points.
	private final DecimalFormat df = createFormatter();

	private DecimalFormat createFormatter() {
		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
		return new DecimalFormat("0.######", symbols);
	}

	private List<String> matrizDosisString = new ArrayList<>();

	// ── Public Utilities ───────────────────────────────────────────────────────

	/**
	 * Extracts Reference and Alternate allele counts from a VCF sample's format field.
	 * Supports GATK (AD), Freebayes (RO/AO), and NGSEP (BSDP).
	 * Dynamically maps A,C,G,T positions for modern NGSEP BSDP tags.
	 * 
	 * @return A double array [refCount, altCount] if counts are found, or null otherwise.
	 */
	public static double[] getRefAltCounts(String[] gData, int adIdx, int bsdpIdx, int roIdx, int aoIdx, String refBase, String altBase) {
		try {
			// 1. NGSEP BSDP Tag (A,C,G,T counts)
			if (bsdpIdx != -1 && gData.length > bsdpIdx && !gData[bsdpIdx].equals(".")) {
				String[] bsdp = gData[bsdpIdx].split(",");
				if (bsdp.length == 4) {
					String firstAlt = altBase.contains(",") ? altBase.split(",")[0] : altBase;
					int refIdx = getBaseIndex(refBase);
					int altIdx = getBaseIndex(firstAlt);
					if (refIdx != -1 && altIdx != -1) {
						return new double[]{Double.parseDouble(bsdp[refIdx]), Double.parseDouble(bsdp[altIdx])};
					}
				} else if (bsdp.length >= 2) {
					// Legacy or simplified BSDP (Alt,Ref)
					return new double[]{Double.parseDouble(bsdp[1]), Double.parseDouble(bsdp[0])};
				}
			} 
			
			// 2. GATK / Standard AD Tag (Ref,Alt1,Alt2...)
			if (adIdx != -1 && gData.length > adIdx && !gData[adIdx].equals(".")) {
				String[] ads = gData[adIdx].split(",");
				if (ads.length >= 2) {
					// Always take the first two: Ref and the primary Alt
					double ref = ads[0].equals(".") ? 0 : Double.parseDouble(ads[0]);
					double alt = ads[1].equals(".") ? 0 : Double.parseDouble(ads[1]);
					return new double[]{ref, alt};
				}
			} 
			
			// 3. Freebayes RO/AO Tags
			if (roIdx != -1 && aoIdx != -1 && gData.length > Math.max(roIdx, aoIdx)) {
				if (!gData[roIdx].equals(".") && !gData[aoIdx].equals(".")) {
					double ref = Double.parseDouble(gData[roIdx]);
					String aoStr = gData[aoIdx];
					double alt = Double.parseDouble(aoStr.contains(",") ? aoStr.split(",")[0] : aoStr);
					return new double[]{ref, alt};
				}
			}
		} catch (Exception ignored) {}
		return null;
	}

	private static int getBaseIndex(String base) {
		if (base == null || base.isEmpty()) return -1;
		switch (base.toUpperCase()) {
			case "A": return 0;
			case "C": return 1;
			case "G": return 2;
			case "T": return 3;
			default: return -1;
		}
	}

	// ── Getters ────────────────────────────────────────────────────────────────

	public int getNumSNPs()      { return numSNPs; }
	public int getNumGenotypes() { return numGenotypes; }

	// ── Public API ─────────────────────────────────────────────────────────────

	/** Backwards compatible method. */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, false, "auto", 0, 5, false);
	}

	/** Backwards compatible method. */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute, boolean storeInMemory) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, "auto", 0, 5, false);
	}

	/** Backwards compatible method. */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, callerType, 0, 5, false);
	}
	
	/** Backwards compatible method. */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType, int minDepth) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, callerType, minDepth, 5, false);
	}
	
	/** Backwards compatible method. */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType, int minDepth, int knnK) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, callerType, minDepth, knnK, false);
	}

	/**
	 * Core function that computes the allele dosage matrix from a VCF file.
	 *
	 * @param vcfFile       Path to the input VCF file containing raw sequence variant calls.
	 * @param ploidy        Ploidy level of the organism. E.g., 10 for sugarcane.
	 * @param impute        Imputation strategy ("bsdp", "mode", "mean", "knn", etc).
	 * @param storeInMemory If true, accumulates the matrix in a List. If false, prints directly to stdout.
	 * @param callerType    The software used to call the variants ("ngsep", "gatk", "freebayes", "auto").
	 * @param minDepth      Statistical Depth Filter. 
	 * @param knnK          Number of neighbors to use if impute strategy contains "knn".
	 * @param adaptiveRounding Whether to use ML 1D Clustering to overcome reference bias.
	 */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType, int minDepth, int knnK, 
			boolean adaptiveRounding) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, callerType, minDepth, knnK, adaptiveRounding, false);
	}

	public List<DosageResult> calculate(String vcfFile, int ploidy) {
		List<DosageResult> results = new ArrayList<>();
		try {
			String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
			numGenotypes = sampleIds.length;
			int n = Math.max(ploidy, 2);
			float[] ploidyLevels = new float[n + 1];
			for (int y = 0; y <= n; y++) ploidyLevels[y] = (1.0f / n) * y;

			Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);
			String lastFormat = "";
			int gtIdx = -1, adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;

			int count = 0;
			for (String[] columns : blockIterator) {
				String chr = columns[0];
				String pos = columns[1];
				String snpId = chr + "_" + pos;

				if (subsetSnps != null && !subsetSnps.contains(snpId)) continue;
				if (count >= maxSnps) break;
				String format = columns.length > 8 ? columns[8] : "";

				if (!format.equals(lastFormat)) {
					lastFormat = format;
					String[] tokens = format.split(":");
					gtIdx = -1; adIdx = -1; roIdx = -1; aoIdx = -1; adpIdx = -1; bsdpIdx = -1;
					for (int i = 0; i < tokens.length; i++) {
						switch (tokens[i]) {
							case "GT":   gtIdx   = i; break;
							case "AD":   adIdx   = i; break;
							case "ADP":  adpIdx  = i; break;
							case "RO":   roIdx   = i; break;
							case "AO":   aoIdx   = i; break;
							case "BSDP": bsdpIdx = i; break;
						}
					}
				}

				DosageResult dr = new DosageResult();
				dr.snpId = chr + "_" + pos;
				dr.dosages = new double[numGenotypes];
				dr.refDepths = new int[numGenotypes];
				dr.altDepths = new int[numGenotypes];
				boolean[] isMissing = new boolean[numGenotypes];
				float[] tmpDosages = new float[numGenotypes];

				extractRawDosagesExtended(columns, numGenotypes, "auto", gtIdx, adIdx, roIdx, aoIdx, adpIdx, bsdpIdx, 
								  0, ploidyLevels, tmpDosages, isMissing, false, true, dr.refDepths, dr.altDepths, null);

				for (int i=0; i<numGenotypes; i++) dr.dosages[i] = tmpDosages[i];
				results.add(dr);
				count++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return results;
	}

	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType, int minDepth, int knnK, 
			boolean adaptiveRounding, boolean rawFrequencies) throws IOException {

		// ── 1. Read sample names from VCF header ───────────────────────────────
		String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
		numGenotypes = sampleIds.length;

		StringBuilder headerBuilder = new StringBuilder("Chr\tpos\t");
		for (String id : sampleIds) {
			headerBuilder.append(id).append("\t");
		}
		emit(headerBuilder.toString(), storeInMemory);

		// ── 2. Pre-compute valid dosage levels for this ploidy ─────────────────
		int n = Math.max(ploidy, 2);
		float[] ploidyLevels = new float[n + 1];
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		// ── Phase 2: K-Nearest Neighbors (KNN) Pass 1 (Optional) ───────────────
		float[][] distanceMatrix = null;
		if (impute.contains("knn")) {
			distanceMatrix = buildDistanceMatrix(vcfFile, callerType, minDepth, ploidyLevels, numGenotypes, adaptiveRounding, 1, "manhattan");
		}

		// ── 3. Stream VCF line by line (low memory footprint) ──────────────────
		System.err.println("[BioJava] Starting allele dosage calculation (Ploidy: " + ploidy + ")...");
		System.err.println("[BioJava] Processing " + numGenotypes + " samples...");
		long startTime = System.currentTimeMillis();
		int snpCount = 0;

		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);

		String lastFormat = "";
		int gtIdx = -1, adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;

		for (String[] columns : blockIterator) {
			snpCount++;
			if (snpCount % 10000 == 0) {
				System.err.print("\r[BioJava] Processed " + snpCount + " variants...");
			}

			String chr    = columns[0];
			String pos    = columns[1];
			String format = columns.length > 8 ? columns[8] : "";

			// ── 4. Parse FORMAT column only when it changes ────────────────────
			if (!format.equals(lastFormat)) {
				lastFormat = format;
				adIdx = -1; roIdx = -1; aoIdx = -1; adpIdx = -1; bsdpIdx = -1;
				String[] tokens = format.split(":");
				for (int i = 0; i < tokens.length; i++) {
					switch (tokens[i]) {
						case "GT":   gtIdx   = i; break;
						case "AD":   adIdx   = i; break;
						case "ADP":  adpIdx  = i; break;
						case "RO":   roIdx   = i; break;
						case "AO":   aoIdx   = i; break;
						case "BSDP": bsdpIdx = i; break;
					}
				}
			}

			// ── 5. Extract raw dosages for all samples in this SNP ─────────────
			int len = Math.min(columns.length - 9, numGenotypes);
			float[] parsedDosages = new float[numGenotypes];
			boolean[] isMissing   = new boolean[numGenotypes];
			
			extractRawDosages(columns, len, callerType, gtIdx, adIdx, roIdx, aoIdx, adpIdx, bsdpIdx, 
							  minDepth, ploidyLevels, parsedDosages, isMissing, adaptiveRounding, rawFrequencies);

			// ── Phase 3: Adaptive Rounding via Clustering (Optional) ───────────
			// If adaptive rounding is enabled, parsedDosages currently holds the purely raw frequencies.
			// We now run K-Means to cluster them and assign them to the biological levels.
			if (adaptiveRounding) {
				assignDosagesViaClustering(parsedDosages, isMissing, ploidyLevels);
			}

			// ── 6. Compute population statistics for non-KNN imputation ────────
			Map<Float, Integer> modeCounter = new HashMap<>();
			for (float nivel : ploidyLevels) {
				modeCounter.put(Float.valueOf(df.format(nivel)), 0);
			}
			float sumDosages     = 0;
			int   numGenotyped   = 0;

			for (int i = 0; i < len; i++) {
				if (!isMissing[i]) {
					modeCounter.merge(parsedDosages[i], 1, Integer::sum);
					sumDosages += parsedDosages[i];
					numGenotyped++;
				}
			}

			// Calculate standard Mean and Mode
			float modeValue = 0;
			if (impute.contains("mode")) {
				int maxCount = -1;
				for (Map.Entry<Float, Integer> e : modeCounter.entrySet()) {
					if (e.getValue() > maxCount) {
						maxCount  = e.getValue();
						modeValue = e.getKey();
					}
				}
			}
			float meanValue = (numGenotyped > 0)
					? Float.valueOf(df.format(sumDosages / numGenotyped))
					: 0;

			// ── 7. Build and emit the output row ──────────────────────────────
			StringBuilder rowBuilder = new StringBuilder();
			rowBuilder.append(chr).append("\t").append(pos).append("\t");

			for (int i = 0; i < len; i++) {
				if (isMissing[i]) {
					if (impute.equals("mode") || impute.equals("bsdp-mode")) {
						rowBuilder.append(modeValue).append("\t");
						
					} else if (impute.equals("mean") || impute.equals("bsdp-mean")) {
						rowBuilder.append(meanValue).append("\t");
						
					} else if (impute.equals("knn") || impute.equals("bsdp-knn")) {
						float knnImputedValue = imputeWithKNN(i, parsedDosages, isMissing, distanceMatrix, knnK);
						knnImputedValue = Float.valueOf(df.format(roundToPloidyLevel(knnImputedValue, ploidyLevels)));
						rowBuilder.append(knnImputedValue).append("\t");
						
					} else { // bsdp (no imputation) → -1
						rowBuilder.append("-1.0\t");
					}
				} else {
					rowBuilder.append(parsedDosages[i]).append("\t");
				}
			}

			emit(rowBuilder.toString(), storeInMemory);
			numSNPs++;
		}
		long endTime = System.currentTimeMillis();
		System.err.println("\n[BioJava] Calculation completed. Total variants exported: " + numSNPs + " | Time: " + (endTime - startTime)/1000.0 + "s");
	}

	// ── Public API: Distance Matrix ────────────────────────────────────────────
	
	public static class DistanceResult {
		public String[] sampleIds;
		public float[][] matrix;
	}

	public DistanceResult computeDistanceMatrix(String vcfFile, int ploidy, String callerType, 
			int minDepth, boolean adaptiveRounding, int threads, String method) throws IOException {
		String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
		int n = Math.max(ploidy, 2);
		float[] ploidyLevels = new float[n + 1];
		for (int y = 0; y <= n; y++) ploidyLevels[y] = (1.0f / n) * y;
		
		float[][] matrix = buildDistanceMatrix(vcfFile, callerType, minDepth, ploidyLevels, sampleIds.length, adaptiveRounding, threads, method);
		
		DistanceResult dr = new DistanceResult();
		dr.sampleIds = sampleIds;
		dr.matrix = matrix;
		return dr;
	}

	/**
	 * Prints the N x N genetic distance matrix between all samples and prints it
	 * as a TSV matrix to the specified output file (or stdout if null).
	 */
	public void computeAndPrintDistanceMatrix(String vcfFile, int ploidy, String callerType, 
			int minDepth, boolean adaptiveRounding, int threads, String outputFile, String method) throws IOException {
		
		java.io.PrintWriter pw;
		if (outputFile != null) {
			pw = new java.io.PrintWriter(new java.io.FileWriter(outputFile));
			System.out.println("[GeneticDistance] Writing distance matrix to: " + outputFile);
		} else {
			pw = new java.io.PrintWriter(System.out);
		}

		String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
		int numG = sampleIds.length;
		
		int n = Math.max(ploidy, 2);
		float[] ploidyLevels = new float[n + 1];
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}
		
		float[][] distanceMatrix = buildDistanceMatrix(vcfFile, callerType, minDepth, ploidyLevels, numG, adaptiveRounding, threads, method);
		
		// Print Header
		pw.print("Sample\t");
		for (int i = 0; i < numG; i++) {
			pw.print(sampleIds[i] + (i == numG - 1 ? "" : "\t"));
		}
		pw.println();
		
		// Print Matrix
		for (int i = 0; i < numG; i++) {
			pw.print(sampleIds[i] + "\t");
			for (int j = 0; j < numG; j++) {
				if (distanceMatrix[i][j] == Float.MAX_VALUE) {
					pw.print("NA" + (j == numG - 1 ? "" : "\t"));
				} else {
					pw.print(df.format(distanceMatrix[i][j]) + (j == numG - 1 ? "" : "\t"));
				}
			}
			pw.println();
		}
		pw.flush();
		if (outputFile != null) pw.close();
	}

	// ── Phase 3 Helpers: Adaptive Rounding via 1D K-Means Clustering ───────────

	/**
	 * Takes raw, unrounded frequencies across the entire population for a single SNP, 
	 * and uses K-Means clustering to assign individuals to true biological states, 
	 * overcoming reference bias.
	 *
	 * Modifies the parsedDosages array in-place with the final rounded values.
	 */
	private void assignDosagesViaClustering(float[] parsedDosages, boolean[] isMissing, float[] ploidyLevels) {
		int k = ploidyLevels.length;
		float[] centroids = Arrays.copyOf(ploidyLevels, k); // Init centers at theoretical levels
		int[] clusterAssignments = new int[parsedDosages.length];
		
		int maxIterations = 20;
		boolean changed;
		
		for (int iter = 0; iter < maxIterations; iter++) {
			changed = false;
			
			// 1. Assign each data point to the closest centroid
			for (int i = 0; i < parsedDosages.length; i++) {
				if (isMissing[i]) continue;
				
				float bestDist = Float.MAX_VALUE;
				int bestCluster = 0;
				for (int c = 0; c < k; c++) {
					float dist = Math.abs(parsedDosages[i] - centroids[c]);
					if (dist < bestDist) {
						bestDist = dist;
						bestCluster = c;
					}
				}
				if (clusterAssignments[i] != bestCluster) {
					clusterAssignments[i] = bestCluster;
					changed = true;
				}
			}
			
			if (!changed && iter > 0) break; // Converged
			
			// 2. Update centroids to the mean of their assigned points
			float[] sum = new float[k];
			int[] count = new int[k];
			
			for (int i = 0; i < parsedDosages.length; i++) {
				if (isMissing[i]) continue;
				int c = clusterAssignments[i];
				sum[c] += parsedDosages[i];
				count[c]++;
			}
			
			for (int c = 0; c < k; c++) {
				if (count[c] > 0) {
					centroids[c] = sum[c] / count[c];
				} else {
					// Empty cluster -> leave it where it was
				}
			}
		}
		
		// 3. Re-assign dosages to the biological (theoretical) level of their cluster
		for (int i = 0; i < parsedDosages.length; i++) {
			if (!isMissing[i]) {
				int c = clusterAssignments[i];
				parsedDosages[i] = Float.valueOf(df.format(ploidyLevels[c]));
			}
		}
	}

	// ── Phase 2 Helpers: KNN Distance Matrix and Imputation ────────────────────
	
	private float[][] buildDistanceMatrix(String vcfFile, String callerType, int minDepth, 
			float[] ploidyLevels, int numGenotypes, boolean adaptiveRounding, int threads, String method) throws IOException {
		
		float[][] distance = new float[numGenotypes][numGenotypes];
		int[][] sharedSnps = new int[numGenotypes][numGenotypes];
		
		// For Nei distance, we need additional matrices to accumulate terms
		final float[][] sumP2 = method.equals("nei") ? new float[numGenotypes][numGenotypes] : null;
		final float[][] sumQ2 = method.equals("nei") ? new float[numGenotypes][numGenotypes] : null;
		final float[][] dotProd = method.equals("nei") ? new float[numGenotypes][numGenotypes] : null;

		// Parallel setup
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		List<Future<?>> futures = new ArrayList<>();
		
		// Per-thread accumulation
		final float[][][] threadDist = (threads > 1) ? new float[threads][numGenotypes][numGenotypes] : null;
		final int[][][] threadShared = (threads > 1) ? new int[threads][numGenotypes][numGenotypes] : null;
		final float[][][] threadSumP2 = (threads > 1 && method.equals("nei")) ? new float[threads][numGenotypes][numGenotypes] : null;
		final float[][][] threadSumQ2 = (threads > 1 && method.equals("nei")) ? new float[threads][numGenotypes][numGenotypes] : null;
		final float[][][] threadDot = (threads > 1 && method.equals("nei")) ? new float[threads][numGenotypes][numGenotypes] : null;
		
		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);
		
		String lastFormat = "";
		int gtIdx = -1, adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;
		int threadCounter = 0;

		for (String[] columns : blockIterator) {
			String formatStr = columns.length > 8 ? columns[8] : "";
			if (!formatStr.equals(lastFormat)) {
				lastFormat = formatStr;
				String[] format = formatStr.split(":");
				gtIdx = -1; adIdx = -1; adpIdx = -1; roIdx = -1; aoIdx = -1; bsdpIdx = -1;
				for (int i = 0; i < format.length; i++) {
					switch (format[i]) {
						case "GT":   gtIdx   = i; break;
						case "AD":   adIdx   = i; break;
						case "ADP":  adpIdx  = i; break;
						case "RO":   roIdx   = i; break;
						case "AO":   aoIdx   = i; break;
						case "BSDP": bsdpIdx = i; break;
					}
				}
			}
			
			final int currentGtIdx = gtIdx, currentAdIdx = adIdx, currentAdpIdx = adpIdx;
			final int currentRoIdx = roIdx, currentAoIdx = aoIdx, currentBsdpIdx = bsdpIdx;
			final String[] currentCols = columns;
			final int currentThreadId = (threads > 1) ? (threadCounter++ % threads) : 0;

			Runnable task = () -> {
				float[] parsedDosages = new float[numGenotypes];
				boolean[] isMissing   = new boolean[numGenotypes];
				int len = Math.min(currentCols.length - 9, numGenotypes);
				
				extractRawDosages(currentCols, len, callerType, currentGtIdx, currentAdIdx, currentRoIdx, currentAoIdx, currentAdpIdx, currentBsdpIdx, 
								  minDepth, ploidyLevels, parsedDosages, isMissing, adaptiveRounding, false);
				
				if (adaptiveRounding) {
					assignDosagesViaClustering(parsedDosages, isMissing, ploidyLevels);
				}

				float[][] targetDist = (threadDist != null) ? threadDist[currentThreadId] : distance;
				int[][] targetShared = (threadShared != null) ? threadShared[currentThreadId] : sharedSnps;
				float[][] targetP2 = (threadSumP2 != null) ? threadSumP2[currentThreadId] : sumP2;
				float[][] targetQ2 = (threadSumQ2 != null) ? threadSumQ2[currentThreadId] : sumQ2;
				float[][] targetDot = (threadDot != null) ? threadDot[currentThreadId] : dotProd;

				for (int i = 0; i < len; i++) {
					if (isMissing[i]) continue;
					for (int j = i + 1; j < len; j++) {
						if (isMissing[j]) continue;
						
						float pi = parsedDosages[i];
						float qi = parsedDosages[j];
						
						switch (method) {
							case "euclidean":
							case "rogers":
								targetDist[i][j] += (pi - qi) * (pi - qi);
								break;
							case "nei":
								targetDot[i][j] += pi * qi;
								targetP2[i][j] += pi * pi;
								targetQ2[i][j] += qi * qi;
								break;
							case "manhattan":
							case "p-distance":
							case "ibs":
							default:
								targetDist[i][j] += Math.abs(pi - qi);
								break;
						}
						targetShared[i][j]++;
					}
				}
			};

			if (threads > 1) {
				futures.add(executor.submit(task));
			} else {
				task.run();
			}
		}
		
		// Wait for completion and merge
		if (threads > 1) {
			for (Future<?> f : futures) {
				try { f.get(); } catch (Exception e) { e.printStackTrace(); }
			}
			executor.shutdown();
			
			// Merge thread-local matrices
			for (int t = 0; t < threads; t++) {
				for (int i = 0; i < numGenotypes; i++) {
					for (int j = i + 1; j < numGenotypes; j++) {
						distance[i][j] += threadDist[t][i][j];
						sharedSnps[i][j] += threadShared[t][i][j];
						if (method.equals("nei")) {
							dotProd[i][j] += threadDot[t][i][j];
							sumP2[i][j] += threadSumP2[t][i][j];
							sumQ2[i][j] += threadSumQ2[t][i][j];
						}
					}
				}
			}
		}

		// Mirro the upper triangle to lower triangle
		for (int i = 0; i < numGenotypes; i++) {
			for (int j = i + 1; j < numGenotypes; j++) {
				distance[j][i] = distance[i][j];
				sharedSnps[j][i] = sharedSnps[i][j];
			}
		}
		
		// Normalize and apply final math transformations
		for (int i = 0; i < numGenotypes; i++) {
			for (int j = 0; j < numGenotypes; j++) {
				if (i == j) {
					distance[i][j] = 0;
					continue;
				}
				
				// Ensure symmetry from the upper triangle
				int row = Math.min(i, j);
				int col = Math.max(i, j);
				
				if (sharedSnps[row][col] > 0) {
					int n = sharedSnps[row][col];
					switch (method) {
						case "euclidean":
							distance[i][j] = (float) Math.sqrt(distance[row][col] / n);
							break;
						case "rogers":
							distance[i][j] = (float) Math.sqrt(distance[row][col] / n); // For bi-allelic it reduces to this
							break;
						case "nei":
							double numerator = dotProd[row][col];
							double denominator = Math.sqrt(sumP2[row][col] * sumQ2[row][col]);
							if (denominator > 0) {
								double identity = numerator / denominator;
								distance[i][j] = (float) -Math.log(Math.max(identity, 1e-10));
							} else {
								distance[i][j] = Float.MAX_VALUE;
							}
							break;
						case "manhattan":
						case "p-distance":
						case "ibs":
						default:
							distance[i][j] = distance[row][col] / n;
							break;
					}
				} else {
					distance[i][j] = Float.MAX_VALUE; // No overlap, infinite distance
				}
			}
		}
		
		return distance;
	}
	
	private static class Neighbor implements Comparable<Neighbor> {
		int index;
		float distance;
		Neighbor(int index, float distance) {
			this.index = index;
			this.distance = distance;
		}
		@Override
		public int compareTo(Neighbor o) {
			return Float.compare(this.distance, o.distance);
		}
	}
	
	private float imputeWithKNN(int targetSampleIdx, float[] parsedDosages, boolean[] isMissing, 
			float[][] distanceMatrix, int k) {
		
		List<Neighbor> validNeighbors = new ArrayList<>();
		
		for (int j = 0; j < parsedDosages.length; j++) {
			if (targetSampleIdx != j && !isMissing[j]) {
				validNeighbors.add(new Neighbor(j, distanceMatrix[targetSampleIdx][j]));
			}
		}
		
		if (validNeighbors.isEmpty()) return 0;
		
		validNeighbors.sort(null);
		
		float sum = 0;
		int count = 0;
		for (int i = 0; i < Math.min(k, validNeighbors.size()); i++) {
			sum += parsedDosages[validNeighbors.get(i).index];
			count++;
		}
		
		return sum / count;
	}

	// ── Core Extractors ────────────────────────────────────────────────────────

	/**
	 * Centralized logic to parse the VCF tags and calculate the statistical dosage.
	 * If adaptiveRounding is TRUE, it leaves the dosages unrounded so the clustering 
	 * algorithm can process them later. Otherwise, it rounds them immediately.
	 */
	private void extractRawDosages(String[] columns, int len, String callerType, 
			int gtIdx, int adIdx, int roIdx, int aoIdx, int adpIdx, int bsdpIdx, int minDepth, 
			float[] ploidyLevels, float[] parsedDosages, boolean[] isMissing, 
			boolean adaptiveRounding, boolean rawFrequencies) {
		extractRawDosagesExtended(columns, len, callerType, gtIdx, adIdx, roIdx, aoIdx, adpIdx, bsdpIdx, minDepth, 
								  ploidyLevels, parsedDosages, isMissing, adaptiveRounding, rawFrequencies, null, null);
	}

	private void extractRawDosagesExtended(String[] columns, int len, String callerType, 
			int gtIdx, int adIdx, int roIdx, int aoIdx, int adpIdx, int bsdpIdx, int minDepth, 
			float[] ploidyLevels, float[] parsedDosages, boolean[] isMissing, 
			boolean adaptiveRounding, boolean rawFrequencies, int[] refDepths, int[] altDepths) {
		
		for (int i = 0; i < len; i++) {
			String   genotypeStr = columns[9 + i];
			String[] gtTokens    = genotypeStr.split(":");

			float   countRef    = 0;
			float   countAlt    = 0;
			boolean foundCounts = false;

			try {
				switch (callerType) {
					case "ngsep":
						if (bsdpIdx != -1 && gtTokens.length > bsdpIdx) {
							String[] bsdp = gtTokens[bsdpIdx].split(",");
							if (bsdp.length == 4 && !bsdp[0].equals(".")) {
								String refBase = columns[3];
								String altBase = columns[4].split(",")[0];
								int refIdx = -1, altIdx = -1;
								if (refBase.equals("A")) refIdx = 0; else if (refBase.equals("C")) refIdx = 1; else if (refBase.equals("G")) refIdx = 2; else if (refBase.equals("T")) refIdx = 3;
								if (altBase.equals("A")) altIdx = 0; else if (altBase.equals("C")) altIdx = 1; else if (altBase.equals("G")) altIdx = 2; else if (altBase.equals("T")) altIdx = 3;
								if (refIdx != -1 && altIdx != -1) {
									countRef = Float.parseFloat(bsdp[refIdx]);
									countAlt = Float.parseFloat(bsdp[altIdx]);
									foundCounts = true;
								}
							} else if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
								countAlt    = Float.parseFloat(bsdp[0]);
								countRef    = Float.parseFloat(bsdp[1]);
								foundCounts = true;
							}
						}
						break;
					case "gatk":
						if (adIdx != -1 && gtTokens.length > adIdx) {
							String[] ads = gtTokens[adIdx].split(",");
							if (ads.length >= 2 && !ads[0].equals(".") && !ads[1].equals(".")) {
								countRef    = Float.parseFloat(ads[0]);
								countAlt    = Float.parseFloat(ads[1]);
								foundCounts = true;
							}
						}
						break;
					case "freebayes":
						if (roIdx != -1 && aoIdx != -1 && gtTokens.length > Math.max(roIdx, aoIdx)) {
							if (!gtTokens[roIdx].equals(".") && !gtTokens[aoIdx].equals(".")) {
								countRef       = Float.parseFloat(gtTokens[roIdx]);
								String aoStr   = gtTokens[aoIdx];
								if (aoStr.contains(",")) aoStr = aoStr.split(",")[0];
								countAlt       = Float.parseFloat(aoStr);
								foundCounts    = true;
							}
						}
						break;
					default: // "auto"
						if (bsdpIdx != -1 && gtTokens.length > bsdpIdx) {
							String[] bsdp = gtTokens[bsdpIdx].split(",");
							if (bsdp.length == 4 && !bsdp[0].equals(".")) {
								String refBase = columns[3];
								String altBase = columns[4].split(",")[0];
								int refIdx = -1, altIdx = -1;
								if (refBase.equals("A")) refIdx = 0; else if (refBase.equals("C")) refIdx = 1; else if (refBase.equals("G")) refIdx = 2; else if (refBase.equals("T")) refIdx = 3;
								if (altBase.equals("A")) altIdx = 0; else if (altBase.equals("C")) altIdx = 1; else if (altBase.equals("G")) altIdx = 2; else if (altBase.equals("T")) altIdx = 3;
								if (refIdx != -1 && altIdx != -1) {
									countRef = Float.parseFloat(bsdp[refIdx]);
									countAlt = Float.parseFloat(bsdp[altIdx]);
									foundCounts = true;
								}
							} else if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
								countAlt = Float.parseFloat(bsdp[0]);
								countRef = Float.parseFloat(bsdp[1]);
								foundCounts = true;
							}
						} else if (adIdx != -1 && gtTokens.length > adIdx) {
							String[] ads = gtTokens[adIdx].split(",");
							if (ads.length >= 2 && !ads[0].equals(".") && !ads[1].equals(".")) {
								countRef = Float.parseFloat(ads[0]);
								countAlt = Float.parseFloat(ads[1]);
								foundCounts = true;
							}
						} else if (adpIdx != -1 && gtTokens.length > adpIdx) {
							String[] adp = gtTokens[adpIdx].split(",");
							if (adp.length >= 2 && !adp[0].equals(".") && !adp[1].equals(".")) {
								countRef = Float.parseFloat(adp[0]);
								countAlt = Float.parseFloat(adp[1]);
								foundCounts = true;
							}
						} else if (roIdx != -1 && aoIdx != -1 && gtTokens.length > Math.max(roIdx, aoIdx)) {
							if (!gtTokens[roIdx].equals(".") && !gtTokens[aoIdx].equals(".")) {
								countRef     = Float.parseFloat(gtTokens[roIdx]);
								String aoStr = gtTokens[aoIdx];
								if (aoStr.contains(",")) aoStr = aoStr.split(",")[0];
								countAlt     = Float.parseFloat(aoStr);
								foundCounts  = true;
							}
						}
				}
				// Capture depths if arrays provided
				if (refDepths != null) refDepths[i] = (int)countRef;
				if (altDepths != null) altDepths[i] = (int)countAlt;

				// Fallback to GT if counts not found
				if (!foundCounts && gtIdx != -1 && gtTokens.length > gtIdx && !gtTokens[gtIdx].startsWith(".")) {
					String[] alleles = gtTokens[gtIdx].split("[/|]");
					for (String a : alleles) {
						if (a.equals("0")) countRef++;
						else if (!a.equals(".")) countAlt++;
					}
					foundCounts = true;
				}
			} catch (NumberFormatException e) {
				// Malformed data → stays missing
			}
			
			// Statistical Depth Filter (Phase 1)
			if (foundCounts && (countRef + countAlt) < minDepth) {
				foundCounts = false;
			}

			isMissing[i] = true;

			if (foundCounts && (countRef + countAlt) > 0) {
				float rawDosage = countAlt / (countRef + countAlt);
				
				if (rawFrequencies) {
					parsedDosages[i] = Float.valueOf(df.format(rawDosage));
				} else if (adaptiveRounding) {
					// Phase 3: Defer rounding. Keep raw frequency.
					parsedDosages[i] = rawDosage;
				} else {
					parsedDosages[i] = Float.valueOf(df.format(roundToPloidyLevel(rawDosage, ploidyLevels)));
				}
				
				isMissing[i] = false;
			}
		}
	}

	// ── Matrix utilities ───────────────────────────────────────────────────────

	/**
	 * Prints the stored dosage matrix to stdout.
	 */
	public void printDosisMatrix() {
		for (String row : matrizDosisString) {
			String[] parts = row.split("\t");
			if (parts.length >= 2) {
				System.out.print(parts[0] + "_" + parts[1] + "\t");
				for (int j = 2; j < parts.length; j++) {
					System.out.print(parts[j] + "\t");
				}
			}
			System.out.println();
		}
	}

	/**
	 * Transposes the matrix so that Samples are rows and SNPs are columns.
	 */
	// ── Internal helpers ───────────────────────────────────────────────────────



	/**
	 * Rounds a raw continuous dosage value to the nearest discrete, valid ploidy level.
	 */
	public float roundToPloidyLevel(float value, float[] ploidyLevels) {
		float bestDist  = Float.MAX_VALUE;
		float rounded   = value;

		for (float level : ploidyLevels) {
			float dist = Math.abs(value - level);
			if (dist < bestDist) {
				bestDist = dist;
				rounded  = level;
			}
		}
		return rounded;
	}

	/** 
	 * Sends a row either to the in-memory list or directly to stdout, facilitating
	 * streaming logic to prevent out-of-memory errors on massive VCF files.
	 */
	private void emit(String row, boolean storeInMemory) {
		if (storeInMemory) {
			matrizDosisString.add(row);
		} else {
			System.out.println(row);
		}
	}

	// ── Per-Individual Ploidy Estimation ───────────────────────────────────────

	/**
	 * Two-pass orchestrator that first estimates the ploidy of each individual
	 * (Pass 1) and then computes the allele dosage matrix using each individual's
	 * specific ploidy for rounding (Pass 2).
	 *
	 * Progress messages are printed to stderr. The ploidy report is printed to
	 * reportStream (or stderr if null) AFTER the dosage matrix is emitted to stdout.
	 *
	 * @param vcfFile            Path to the input VCF file.
	 * @param candidatePloidies  Array of ploidy values to evaluate (e.g. {2,4,6,8,10,12}).
	 * @param impute             Imputation strategy.
	 * @param storeInMemory      If true, stores rows in memory instead of printing to stdout.
	 * @param callerType         Variant caller type.
	 * @param minDepth           Minimum read depth filter.
	 * @param knnK               KNN neighbors (only relevant if impute contains "knn").
	 * @param adaptiveRounding   Enable K-Means adaptive rounding.
	 * @param rawFrequencies     Export raw frequencies instead of rounded dosages.
	 * @param reportStream       Where to print the ploidy report (null → stderr).
	 */


}
