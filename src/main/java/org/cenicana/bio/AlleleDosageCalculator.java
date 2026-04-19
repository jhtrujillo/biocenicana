package org.cenicana.bio;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlleleDosageCalculator {

	private int numSNPs = 0;
	private int numGenotypes = 0;
	private final DecimalFormat df = new DecimalFormat("#.0");

	private String[][] matrizTranspuestaDosis;
	private List<String> matrizDosisString = new ArrayList<>();

	// ── Getters ────────────────────────────────────────────────────────────────

	public int getNumSNPs()      { return numSNPs; }
	public int getNumGenotypes() { return numGenotypes; }

	// ── Public API ─────────────────────────────────────────────────────────────

	public void computeAlleleDosage(String vcfFile, int ploidy, String impute) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, false, "auto");
	}

	public void computeAlleleDosage(String vcfFile, int ploidy, String impute, boolean storeInMemory) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, "auto");
	}

	/**
	 * Computes an allele dosage matrix from a VCF file.
	 *
	 * @param vcfFile       Path to the VCF file.
	 * @param ploidy        Ploidy level of the organism (e.g. 10 for sugarcane).
	 * @param impute        Imputation strategy: "bsdp" | "mode" | "mean" | "bsdp-mode" | "bsdp-mean".
	 * @param storeInMemory If true, stores the matrix in memory instead of printing to stdout.
	 * @param callerType    Variant caller that generated the VCF: "ngsep" | "gatk" | "freebayes" | "auto".
	 */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute,
			boolean storeInMemory, String callerType) throws IOException {

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

		// ── 3. Stream VCF line by line (low memory) ────────────────────────────
		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);

		// Cache FORMAT parsing: most VCFs use the same FORMAT for every SNP.
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

			// ── 5. Initialise per-SNP accumulators (reuse ploidyLevels) ────────
			Map<Float, Integer> modeCounter = new HashMap<>();
			for (float nivel : ploidyLevels) {
				modeCounter.put(Float.valueOf(df.format(nivel)), 0);
			}
			float sumDosages     = 0;
			int   numGenotyped   = 0;

			int len = Math.min(columns.length - 9, sampleIds.length);
			float[] parsedDosages = new float[sampleIds.length];
			boolean[] isMissing   = new boolean[sampleIds.length];

			// ── 6. Compute dosage per sample ───────────────────────────────────
			for (int i = 0; i < len; i++) {
				String   genotypeStr = columns[9 + i];
				String[] gtTokens    = genotypeStr.split(":");

				float   countRef    = 0;
				float   countAlt    = 0;
				boolean foundCounts = false;

				try {
					switch (callerType) {
						case "ngsep":
							// NGSEP BSDP format: [countAlt, countRef, ...]
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
							// GATK AD format: [countRef, countAlt]
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
							// FreeBayes: RO (ref) + AO (alt)
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

						default: // "auto" — try in priority order: BSDP → AD → ADP → RO+AO
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

				isMissing[i] = true;

				if (foundCounts && (countRef + countAlt) > 0) {
					float rawDosage = countRef / (countRef + countAlt);
					float dosage    = Float.valueOf(df.format(roundToPloidyLevel(rawDosage, ploidyLevels)));

					modeCounter.merge(dosage, 1, Integer::sum);
					parsedDosages[i] = dosage;
					isMissing[i]     = false;
					sumDosages       += dosage;
					numGenotyped++;
				}
			}

			// ── 7. Compute imputation values ───────────────────────────────────
			float modeValue = 0;
			if (impute.equals("mode") || impute.equals("bsdp-mode")) {
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

			// ── 8. Build and emit the output row ──────────────────────────────
			StringBuilder rowBuilder = new StringBuilder();
			rowBuilder.append(chr).append("\t").append(pos).append("\t");

			for (int i = 0; i < len; i++) {
				if (isMissing[i]) {
					if (impute.equals("mode") || impute.equals("bsdp-mode")) {
						rowBuilder.append(modeValue).append("\t");
					} else if (impute.equals("mean") || impute.equals("bsdp-mean")) {
						rowBuilder.append(meanValue).append("\t");
					} else { // bsdp → -1
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

	// ── Matrix utilities ───────────────────────────────────────────────────────

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
	 * Rounds a dosage value to the nearest valid ploidy level.
	 *
	 * @param value       Raw dosage ratio (0.0 – 1.0).
	 * @param ploidyLevels Array of valid dosage levels for the organism's ploidy.
	 * @return The closest value in ploidyLevels.
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

	/** Sends a row either to the in-memory list or directly to stdout. */
	private void emit(String row, boolean storeInMemory) {
		if (storeInMemory) {
			matrizDosisString.add(row);
		} else {
			System.out.println(row);
		}
	}
}