package org.cenicana.bio;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
	
	// Formatter to round decimal representations to a single decimal point.
	private final DecimalFormat df = new DecimalFormat("#.0");

	private String[][] matrizTranspuestaDosis;
	private List<String> matrizDosisString = new ArrayList<>();

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
			distanceMatrix = buildDistanceMatrix(vcfFile, callerType, minDepth, ploidyLevels, numGenotypes, adaptiveRounding);
		}

		// ── 3. Stream VCF line by line (low memory footprint) ──────────────────
		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);

		String lastFormat = "";
		int adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;

		for (String[] columns : blockIterator) {

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
			
			extractRawDosages(columns, len, callerType, adIdx, roIdx, aoIdx, adpIdx, bsdpIdx, 
							  minDepth, ploidyLevels, parsedDosages, isMissing, adaptiveRounding);

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
	}

	// ── Public API: Distance Matrix ────────────────────────────────────────────

	/**
	 * Computes the N x N genetic distance matrix between all samples and prints it
	 * as a TSV matrix to stdout.
	 */
	public void computeAndPrintDistanceMatrix(String vcfFile, int ploidy, String callerType, 
			int minDepth, boolean adaptiveRounding) throws IOException {
		
		String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
		int numG = sampleIds.length;
		
		int n = Math.max(ploidy, 2);
		float[] ploidyLevels = new float[n + 1];
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}
		
		float[][] distanceMatrix = buildDistanceMatrix(vcfFile, callerType, minDepth, ploidyLevels, numG, adaptiveRounding);
		
		// Print Header
		System.out.print("Sample\t");
		for (int i = 0; i < numG; i++) {
			System.out.print(sampleIds[i] + (i == numG - 1 ? "" : "\t"));
		}
		System.out.println();
		
		// Print Matrix
		for (int i = 0; i < numG; i++) {
			System.out.print(sampleIds[i] + "\t");
			for (int j = 0; j < numG; j++) {
				if (distanceMatrix[i][j] == Float.MAX_VALUE) {
					System.out.print("NA" + (j == numG - 1 ? "" : "\t"));
				} else {
					System.out.print(df.format(distanceMatrix[i][j]) + (j == numG - 1 ? "" : "\t"));
				}
			}
			System.out.println();
		}
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
			float[] ploidyLevels, int numGenotypes, boolean adaptiveRounding) throws IOException {
		
		float[][] distance = new float[numGenotypes][numGenotypes];
		int[][] sharedSnps = new int[numGenotypes][numGenotypes];
		
		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);
		
		String lastFormat = "";
		int adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;
		
		float[] parsedDosages = new float[numGenotypes];
		boolean[] isMissing   = new boolean[numGenotypes];
		
		for (String[] columns : blockIterator) {
			String format = columns.length > 8 ? columns[8] : "";
			if (!format.equals(lastFormat)) {
				lastFormat = format;
				adIdx = -1; roIdx = -1; aoIdx = -1; adpIdx = -1; bsdpIdx = -1;
				String[] tokens = format.split(":");
				for (int i = 0; i < tokens.length; i++) {
					switch (tokens[i]) {
						case "AD":   adIdx   = i; break;
						case "ADP":  adpIdx  = i; break;
						case "RO":   roIdx   = i; break;
						case "AO":   aoIdx   = i; break;
						case "BSDP": bsdpIdx = i; break;
					}
				}
			}
			
			int len = Math.min(columns.length - 9, numGenotypes);
			extractRawDosages(columns, len, callerType, adIdx, roIdx, aoIdx, adpIdx, bsdpIdx, 
							  minDepth, ploidyLevels, parsedDosages, isMissing, adaptiveRounding);
			
			if (adaptiveRounding) {
				assignDosagesViaClustering(parsedDosages, isMissing, ploidyLevels);
			}
			
			// Accumulate distance
			for (int i = 0; i < len; i++) {
				if (isMissing[i]) continue;
				for (int j = i + 1; j < len; j++) {
					if (isMissing[j]) continue;
					distance[i][j] += Math.abs(parsedDosages[i] - parsedDosages[j]);
					distance[j][i] = distance[i][j];
					sharedSnps[i][j]++;
					sharedSnps[j][i]++;
				}
			}
		}
		
		// Normalize by number of shared SNPs to prevent bias from missing data
		for (int i = 0; i < numGenotypes; i++) {
			for (int j = 0; j < numGenotypes; j++) {
				if (i != j && sharedSnps[i][j] > 0) {
					distance[i][j] /= sharedSnps[i][j];
				} else if (i != j) {
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
			int adIdx, int roIdx, int aoIdx, int adpIdx, int bsdpIdx, int minDepth, 
			float[] ploidyLevels, float[] parsedDosages, boolean[] isMissing, boolean adaptiveRounding) {
		
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
							if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
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
							if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
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
						break;
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
				float rawDosage = countRef / (countRef + countAlt);
				
				if (adaptiveRounding) {
					// Phase 3: Defer rounding. Keep raw frequency.
					parsedDosages[i] = rawDosage;
				} else {
					// Legacy: Strict mathematical rounding
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
	public void TransposeDosisMatrix() {
		int numRows = numGenotypes + 2;
		int numCol  = numSNPs + 1;
		matrizTranspuestaDosis = new String[numRows][numCol];

		for (int i = 0; i < numCol; i++) {
			String[] tmp = matrizDosisString.get(i).split("\t");
			for (int j = 0; j < tmp.length; j++) {
				matrizTranspuestaDosis[j][i] = tmp[j];
			}
		}
	}

	/**
	 * Prints the transposed matrix to stdout.
	 */
	public void printTransposeDosisMatrix() {
		for (String[] row : matrizTranspuestaDosis) {
			for (String cell : row) {
				System.out.print(cell + "\t");
			}
			System.out.println();
		}
	}

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

