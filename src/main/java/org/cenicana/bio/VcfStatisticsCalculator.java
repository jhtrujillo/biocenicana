package org.cenicana.bio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import org.cenicana.bio.io.VcfFastReader;

public class VcfStatisticsCalculator {

	// ── Global variant counts ────────────────────────────────────────────────
	public int numSnps       = 0;
	public int numIndels     = 0;
	public int numTransitions   = 0;
	public int numTransversions = 0;

	// ── Per-sample metrics ───────────────────────────────────────────────────
	public String[] sampleNames;
	public int[]    sampleMissingCount;   // # missing genotypes per sample
	public int[]    sampleGenotypedCount; // # valid (non-missing) genotypes per sample
	public int[]    sampleHetCount;       // # heterozygous genotypes per sample
	public long[]   sampleTotalDepth;     // sum of DP values per sample
	public int[]    sampleDepthCount;     // # genotypes that had a DP value

	// Derived per-sample statistics (calculated after the loop)
	public double[] sampleObsHet;   // Observed Heterozygosity (OH) = hetCount/genotypedCount
	public double[] sampleFstat;    // F-statistic = 1 - OH/mean_EH

	// ── Global EH accumulator (to compute mean EH across all sites) ──────────
	private double ehSum   = 0;
	private int    ehCount = 0;
	public  double meanEH  = 0;  // Mean Expected Heterozygosity across all SNPs

	// ── Variant density per chromosome ───────────────────────────────────────
	// LinkedHashMap preserves insertion order (chromosome order)
	public Map<String, Integer> variantsPerChromosome = new LinkedHashMap<>();

	// ── Histograms ───────────────────────────────────────────────────────────
	// Site missingness bins: [0-10%, 10-20%, ... 90-100%] = 10 bins
	public int[] siteMissingnessHistogram = new int[10];

	// MAF bins: [0-0.05, 0.05-0.10, ... 0.45-0.50] = 10 bins
	public int[] mafHistogram = new int[10];

	// EH per-site histogram: [0-0.1, 0.1-0.2, ... 0.9-1.0] = 10 bins
	public int[] ehHistogram = new int[10];

	// ── Main calculation method ───────────────────────────────────────────────
	public void calculate(String vcfFile) throws IOException {
		sampleNames = VcfFastReader.getSampleIds(vcfFile);
		int numSamples = sampleNames.length;

		sampleMissingCount   = new int[numSamples];
		sampleGenotypedCount = new int[numSamples];
		sampleHetCount       = new int[numSamples];
		sampleTotalDepth     = new long[numSamples];
		sampleDepthCount     = new int[numSamples];

		try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) continue;

				String[] cols = line.split("\t");
				if (cols.length < 9 + numSamples) continue;

				String chrom = cols[0];
				String ref   = cols[3];
				String alt   = cols[4];

				// ── Variant type ──────────────────────────────────────────────
				boolean isBiallelicSnp = ref.length() == 1 && alt.length() == 1 && !alt.equals(".");
				if (isBiallelicSnp) {
					numSnps++;
					if (isTransition(ref, alt)) numTransitions++;
					else                         numTransversions++;
				} else {
					numIndels++;
				}

				// ── Density per chromosome ────────────────────────────────────
				variantsPerChromosome.merge(chrom, 1, Integer::sum);

				// ── Format field parsing ──────────────────────────────────────
				String[] formatFields = cols[8].split(":");
				int depthIndex = -1;
				for (int i = 0; i < formatFields.length; i++) {
					if (formatFields[i].equals("BSDP") || formatFields[i].equals("DP")) {
						depthIndex = i;
						break;
					}
				}

				// ── Allele counters for this site ─────────────────────────────
				int siteMissing = 0;
				int totalAlleles = 0;
				int refAlleles   = 0;
				int altAlleles   = 0;

				// ── Genotype loop ─────────────────────────────────────────────
				for (int i = 0; i < numSamples; i++) {
					String gData = cols[9 + i];
					String gt    = gData.split(":")[0];

					if (gt.equals(".") || gt.startsWith("./.") || gt.startsWith("./")) {
						sampleMissingCount[i]++;
						siteMissing++;
						continue;
					}

					sampleGenotypedCount[i]++;

					// Parse allele calls (handles / and | separators, and polyploids)
					String[] alleles = gt.split("[/|]");
					boolean hasRef = false, hasAlt = false;
					for (String a : alleles) {
						if (a.equals(".")) continue;
						totalAlleles++;
						if (a.equals("0")) { refAlleles++; hasRef = true; }
						else               { altAlleles++; hasAlt = true; }
					}
					// Heterozygous = has both ref and alt alleles
					if (hasRef && hasAlt) sampleHetCount[i]++;

					// ── Depth ─────────────────────────────────────────────────
					if (depthIndex != -1) {
						String[] sFields = gData.split(":");
						if (sFields.length > depthIndex && !sFields[depthIndex].equals(".")) {
							try {
								// BSDP is comma-separated (A,C,G,T), sum all
								String depthRaw = sFields[depthIndex];
								int depth = 0;
								if (depthRaw.contains(",")) {
									for (String d : depthRaw.split(",")) depth += Integer.parseInt(d.trim());
								} else {
									depth = Integer.parseInt(depthRaw);
								}
								sampleTotalDepth[i] += depth;
								sampleDepthCount[i]++;
							} catch (NumberFormatException ignored) {}
						}
					}
				}

				// ── Site-level histograms ─────────────────────────────────────
				double siteMissFreq = (double) siteMissing / numSamples;
				int binMissing = Math.min((int)(siteMissFreq * 10), 9);
				siteMissingnessHistogram[binMissing]++;

				if (totalAlleles > 0) {
					double altFreq = (double) altAlleles / totalAlleles;
					double maf     = Math.min(altFreq, 1.0 - altFreq);

					// MAF histogram
					int binMaf = Math.min((int)(maf / 0.05), 9);
					mafHistogram[binMaf]++;

					// Expected Heterozygosity EH = 1 - (p^2 + q^2) = 2pq
					double p = (double) refAlleles / totalAlleles;
					double q = 1.0 - p;
					double eh = 2.0 * p * q;
					ehSum += eh;
					ehCount++;

					// EH histogram
					int binEh = Math.min((int)(eh * 10), 9);
					ehHistogram[binEh]++;
				}
			}
		}

		// ── Derived calculations (post-loop) ──────────────────────────────────
		meanEH = ehCount > 0 ? ehSum / ehCount : 0.0;

		sampleObsHet = new double[sampleNames.length];
		sampleFstat  = new double[sampleNames.length];
		for (int i = 0; i < sampleNames.length; i++) {
			double oh = sampleGenotypedCount[i] > 0
				? (double) sampleHetCount[i] / sampleGenotypedCount[i]
				: 0.0;
			sampleObsHet[i] = oh;
			// F = 1 - OH/EH; if meanEH == 0, define F = 0
			sampleFstat[i]  = meanEH > 0 ? 1.0 - (oh / meanEH) : 0.0;
		}
	}

	// ── Helpers ───────────────────────────────────────────────────────────────
	private boolean isTransition(String ref, String alt) {
		ref = ref.toUpperCase();
		alt = alt.toUpperCase();
		return (ref.equals("A") && alt.equals("G")) ||
			   (ref.equals("G") && alt.equals("A")) ||
			   (ref.equals("C") && alt.equals("T")) ||
			   (ref.equals("T") && alt.equals("C"));
	}
}


	public void calculate(String vcfFile) throws IOException {
		sampleNames = VcfFastReader.getSampleIds(vcfFile);
		int numSamples = sampleNames.length;

		sampleMissingCount = new int[numSamples];
		sampleTotalDepth = new long[numSamples];
		sampleDepthCount = new int[numSamples];

		try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					continue;
				}

				String[] cols = line.split("\t");
				if (cols.length < 9 + numSamples) {
					continue; // Invalid row
				}

				String ref = cols[3];
				String alt = cols[4];
				
				// Identify SNP vs InDel
				if (ref.length() == 1 && alt.length() == 1 && !alt.equals(".")) {
					numSnps++;
					// Ts/Tv calculation
					if (isTransition(ref, alt)) {
						numTransitions++;
					} else {
						numTransversions++;
					}
				} else {
					numIndels++;
				}

				// Find DP/BSDP index in format field
				String format = cols[8];
				String[] formatFields = format.split(":");
				int depthIndex = -1;
				for (int i = 0; i < formatFields.length; i++) {
					if (formatFields[i].equals("BSDP") || formatFields[i].equals("DP")) {
						depthIndex = i;
						break;
					}
				}

				int siteMissing = 0;
				int totalAlleles = 0;
				int altAlleles = 0;

				// Parse genotypes
				for (int i = 0; i < numSamples; i++) {
					String genotypeData = cols[9 + i];
					
					if (genotypeData.startsWith("./.") || genotypeData.startsWith(".")) {
						sampleMissingCount[i]++;
						siteMissing++;
					} else {
						// Calculate MAF (very basic count based on 0/1 splits)
						String[] gtParts = genotypeData.split(":")[0].split("[/|]");
						for (String allele : gtParts) {
							if (!allele.equals(".")) {
								totalAlleles++;
								if (!allele.equals("0")) {
									altAlleles++;
								}
							}
						}

						// Calculate Depth
						if (depthIndex != -1) {
							String[] sampleFields = genotypeData.split(":");
							if (sampleFields.length > depthIndex && !sampleFields[depthIndex].equals(".")) {
								try {
									int depth = Integer.parseInt(sampleFields[depthIndex]);
									sampleTotalDepth[i] += depth;
									sampleDepthCount[i]++;
								} catch (NumberFormatException e) {
									// Ignore parsing error for depth
								}
							}
						}
					}
				}

				// Update Histograms
				double siteMissFreq = (double) siteMissing / numSamples;
				int binMissing = (int) (siteMissFreq * 10);
				if (binMissing >= 10) binMissing = 9;
				siteMissingnessHistogram[binMissing]++;

				if (totalAlleles > 0) {
					double altFreq = (double) altAlleles / totalAlleles;
					double maf = Math.min(altFreq, 1.0 - altFreq);
					int binMaf = (int) (maf * 20); // 0.0 to 0.5 -> 20 bins (size 0.025), wait 0 to 0.5 needs 10 bins of 0.05
					binMaf = (int) (maf / 0.05);
					if (binMaf >= 10) binMaf = 9;
					mafHistogram[binMaf]++;
				}
			}
		}
	}

	private boolean isTransition(String ref, String alt) {
		ref = ref.toUpperCase();
		alt = alt.toUpperCase();
		return (ref.equals("A") && alt.equals("G")) ||
				(ref.equals("G") && alt.equals("A")) ||
				(ref.equals("C") && alt.equals("T")) ||
				(ref.equals("T") && alt.equals("C"));
	}
}
