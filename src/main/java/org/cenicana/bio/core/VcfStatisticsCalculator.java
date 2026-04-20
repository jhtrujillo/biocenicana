package org.cenicana.bio.core;
import org.cenicana.bio.utils.FileUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.cenicana.bio.io.VcfFastReader;

/**
 * Streaming VCF statistics calculator.
 * Level 1: Ts/Tv, MAF, Missingness, Depth, OH, F-statistic, EH, Density per chromosome.
 * Level 2: HWE chi-square & Fisher Exact Test, AN distribution, Tajima's D, pairwise Fst.
 */
public class VcfStatisticsCalculator {

	// ── Global variant counts ────────────────────────────────────────────────
	public int numSnps          = 0;
	public int numIndels        = 0;
	public int numTransitions   = 0;
	public int numTransversions = 0;

	// ── Per-sample metrics ───────────────────────────────────────────────────
	public String[] sampleNames;
	public int[]    sampleMissingCount;    // # missing genotypes per sample
	public int[]    sampleGenotypedCount;  // # valid (non-missing) genotypes per sample
	public int[]    sampleHetCount;        // # heterozygous genotypes per sample
	public int[]    sampleHomoRefCount;    // # homozygous reference genotypes
	public int[]    sampleHomoAltCount;    // # homozygous alternative genotypes
	public int[]    sampleNonRefCount;     // # any genotype with >=1 alt allele
	public int[]    sampleTsCount;         // # transitions in this sample
	public int[]    sampleTvCount;         // # transversions in this sample
	public int[]    sampleRareAlleleCount; // # rare alleles (MAF < 0.05) in this sample
	public long[]   sampleTotalDepth;      // sum of DP values per sample
	public int[]    sampleDepthCount;      // # genotypes that had a DP value

	// Derived per-sample statistics (post-loop)
	public double[] sampleObsHet;    // OH = hetCount / genotypedCount
	public double[] sampleFstat;     // F = 1 - OH / meanEH

	// ── Global EH accumulator ────────────────────────────────────────────────
	private double ehSum   = 0;
	private int    ehCount = 0;
	public  double meanEH  = 0;

	// ── Density per chromosome ────────────────────────────────────────────────
	public Map<String, Integer> variantsPerChromosome = new LinkedHashMap<>();
	// Bin size for genomic density (1Mb)
	public static final int DENSITY_BIN_SIZE = 1000000;
	// Map<Chrom, Map<BinIdx, Count>>
	public Map<String, Map<Integer, Integer>> binnedDensity = new LinkedHashMap<>();

	// ── Histograms ────────────────────────────────────────────────────────────
	public int[] siteMissingnessHistogram = new int[10]; // 0-10% ... 90-100%
	public int[] mafHistogram             = new int[50]; // 0-0.01 ... 0.49-0.50 (50 bins)
	public int[] ehHistogram              = new int[10]; // 0.0-0.1 ... 0.9-1.0

	// ════════════════════════════════════════════════════════════════════════
	// LEVEL 2 fields
	// ════════════════════════════════════════════════════════════════════════

	// ── HWE (biallelic diploid SNPs only) ────────────────────────────────────
	// bins: [0-2), [2-3.84), [3.84-7), [7-10), [10+)
	public int[]  hweChiSqHistogram  = new int[5];
	public int    numHweViolations   = 0;  // sites with p < 0.05 (from Fisher exact test)
	public int    numHweTested       = 0;  // biallelic diploid SNPs with enough data
	public double meanChiSq          = 0;
	private double chiSqSum          = 0;

	// ── AN – distinct alleles per site ───────────────────────────────────────
	public int numMonomorphic  = 0;  // no alt allele observed
	public int numBiallelic    = 0;  // exactly one alt allele
	public int numMultiallelic = 0;  // 2+ alt alleles in ALT column

	// ── Tajima's D (genome-wide, diploid assumption) ─────────────────────────
	public double tajimaD       = Double.NaN;
	public int    numSegSites   = 0;   // S: segregating sites
	public double piHat         = 0;   // π: average pairwise differences per site
	public double thetaW        = 0;   // Watterson's θ_W = S / a1
	private double piPerSiteSum = 0;   // accumulated sum of per-site π_i
	private long   totalAllelesForTajima = 0; // for computing effective n
	private int    tajimaSiteCount       = 0; // sites used for Tajima's

	// ── Fst between populations (optional) ───────────────────────────────────
	private Map<String, Integer> samplePopIndex = null; // sample name → pop index
	public  String[]             populationNames = null;
	public  double[][]           pairwiseFst    = null; // [pop_i][pop_j] mean Fst

	/** Call this BEFORE calculate() to enable population-level Fst analysis. */
	public void loadPopulationMap(String popFile) throws IOException {
		Map<String, String> sampleToPop = new HashMap<>();
		List<String> popOrder = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(popFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) continue;
				String[] parts = line.split("\t");
				if (parts.length < 2) continue;
				String sample = parts[0].trim();
				String pop    = parts[1].trim();
				sampleToPop.put(sample, pop);
				if (!popOrder.contains(pop)) popOrder.add(pop);
			}
		}
		populationNames  = popOrder.toArray(new String[0]);
		samplePopIndex   = new HashMap<>();
		for (int i = 0; i < populationNames.length; i++) {
			final int idx = i;
			sampleToPop.forEach((s, p) -> {
				if (p.equals(populationNames[idx])) samplePopIndex.put(s, idx);
			});
		}
	}

	// ── Main calculation method ───────────────────────────────────────────────
	public void calculate(String vcfFile) throws IOException {
		sampleNames = VcfFastReader.getSampleIds(vcfFile);
		int numSamples = sampleNames.length;

		sampleMissingCount    = new int[numSamples];
		sampleGenotypedCount  = new int[numSamples];
		sampleHetCount        = new int[numSamples];
		sampleHomoRefCount    = new int[numSamples];
		sampleHomoAltCount    = new int[numSamples];
		sampleNonRefCount     = new int[numSamples];
		sampleTsCount         = new int[numSamples];
		sampleTvCount         = new int[numSamples];
		sampleRareAlleleCount = new int[numSamples];
		sampleTotalDepth      = new long[numSamples];
		sampleDepthCount      = new int[numSamples];

		// Build sample→pop index array
		int[] samplePopIdx = null;
		int   numPops      = 0;
		double[][][] fstAccum = null;
		if (samplePopIndex != null && populationNames != null) {
			numPops      = populationNames.length;
			samplePopIdx = new int[numSamples];
			fstAccum     = new double[numPops][numPops][3];
			for (int i = 0; i < numSamples; i++) {
				Integer pidx = samplePopIndex.get(sampleNames[i]);
				samplePopIdx[i] = (pidx != null) ? pidx : -1;
			}
		}

		try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
			String line;
			int[][] parsedGenotypes = new int[numSamples][]; // Temp array for 2-pass per site

			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) continue;

				String[] cols = line.split("\t");
				if (cols.length < 9 + numSamples) continue;

				String chrom = cols[0];
				int    pos   = Integer.parseInt(cols[1]);
				String ref   = cols[3];
				String alt   = cols[4];
				
				// Update counters
				variantsPerChromosome.merge(chrom, 1, Integer::sum);
				int binIdx = pos / DENSITY_BIN_SIZE;
				binnedDensity.computeIfAbsent(chrom, k -> new HashMap<>()).merge(binIdx, 1, Integer::sum);
				String[] altAllelesList = alt.equals(".") ? new String[0] : alt.split(",");
				int numAllelesAtSite = 1 + altAllelesList.length;

				if (altAllelesList.length == 0) numMonomorphic++;
				else if (altAllelesList.length == 1) numBiallelic++;
				else numMultiallelic++;

				boolean isBiallelicSnp = ref.length() == 1
					&& altAllelesList.length == 1
					&& altAllelesList[0].length() == 1
					&& !alt.equals(".");
				boolean isTransition = isBiallelicSnp && isTransition(ref, altAllelesList[0]);

				if (isBiallelicSnp) {
					numSnps++;
					if (isTransition) numTransitions++;
					else              numTransversions++;
				} else {
					numIndels++;
				}

				String[] fmt  = cols[8].split(":");
				int depthIdx  = -1;
				for (int i = 0; i < fmt.length; i++) {
					if (fmt[i].equals("BSDP") || fmt[i].equals("DP")) { depthIdx = i; break; }
				}

				// PASS 1: Parse alleles & counts
				int siteMissing = 0;
				int totalGenotypedAlleles = 0;
				int[] siteAlleleCounts = new int[numAllelesAtSite];
				int n00 = 0, n01 = 0, n11 = 0; // HWE stats
				int[][] popAlleles = (samplePopIdx != null) ? new int[numPops][2] : null;

				for (int i = 0; i < numSamples; i++) {
					parsedGenotypes[i] = null;
					String gData = cols[9 + i];
					String gt    = gData.split(":")[0];

					if (gt.equals(".") || gt.startsWith("./.") || gt.startsWith("./")) {
						siteMissing++;
					} else {
						String[] alleles = gt.split("[/|]");
						int[] numericAlleles = new int[alleles.length];
						boolean valid = true;
						for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
							if (alleles[aIdx].equals(".")) {
								valid = false;
								break;
							}
							try {
								numericAlleles[aIdx] = Integer.parseInt(alleles[aIdx]);
								if (numericAlleles[aIdx] >= numAllelesAtSite) valid = false; // invalid VCF GT
							} catch (NumberFormatException e) {
								valid = false;
							}
						}

						if (valid) {
							parsedGenotypes[i] = numericAlleles;
							sampleGenotypedCount[i]++;
							for (int a : numericAlleles) {
								siteAlleleCounts[a]++;
								totalGenotypedAlleles++;
							}

							// HWE counts (diploid biallelic)
							if (isBiallelicSnp && numericAlleles.length == 2) {
								if (numericAlleles[0] == 0 && numericAlleles[1] == 0) n00++;
								else if (numericAlleles[0] == 1 && numericAlleles[1] == 1) n11++;
								else n01++;
							}

							// Pop allele counts (biallelic only)
							if (popAlleles != null && samplePopIdx[i] >= 0 && isBiallelicSnp) {
								for (int a : numericAlleles) {
									popAlleles[samplePopIdx[i]][a]++;
								}
							}
						} else {
							siteMissing++; // fallback to missing
						}
					}

					if (depthIdx != -1) {
						String[] sf = gData.split(":");
						if (sf.length > depthIdx && !sf[depthIdx].equals(".")) {
							try {
								String dr = sf[depthIdx];
								int depth = 0;
								if (dr.contains(",")) {
									for (String d : dr.split(",")) depth += Integer.parseInt(d.trim());
								} else {
									depth = Integer.parseInt(dr);
								}
								sampleTotalDepth[i] += depth;
								sampleDepthCount[i]++;
							} catch (NumberFormatException ignored) {}
						}
					}
				}

				// Site missingness
				double missFreq = (double) siteMissing / numSamples;
				siteMissingnessHistogram[Math.min((int)(missFreq * 10), 9)]++;

				if (totalGenotypedAlleles > 0) {
					// Expected Heterozygosity: 1 - sum(pi^2) (Handles multiallelic accurately like NGSEP)
					double sumSqFreqs = 0;
					int minAc = Integer.MAX_VALUE;
					int minAcIdx = -1;
					for (int j = 0; j < siteAlleleCounts.length; j++) {
						int c = siteAlleleCounts[j];
						if (c > 0) {
							double freq = (double) c / totalGenotypedAlleles;
							sumSqFreqs += freq * freq;
							if (c < minAc) {
								minAc = c;
								minAcIdx = j;
							}
						}
					}
					double eh = 1.0 - sumSqFreqs;
					ehSum += eh;
					ehCount++;
					ehHistogram[Math.min((int)(eh * 10), 9)]++;

					// MAF
					double maf = (double) minAc / totalGenotypedAlleles;
					mafHistogram[Math.min((int)(maf * 100), 49)]++; // 50 bins

					// PASS 2: Per-sample metrics using MAF & exact genotypes
					for (int i = 0; i < numSamples; i++) {
						int[] gt = parsedGenotypes[i];
						if (gt == null) {
							sampleMissingCount[i]++;
							continue;
						}

						boolean isHet = false;
						boolean hasRef = false;
						boolean hasAlt = false;
						boolean hasRare = false;

						for (int a : gt) {
							if (a == 0) hasRef = true;
							else hasAlt = true;
							if (a == minAcIdx && maf < 0.05) hasRare = true;
						}

						// Check if all alleles are the same
						for (int j = 1; j < gt.length; j++) {
							if (gt[j] != gt[0]) {
								isHet = true;
								break;
							}
						}

						if (isHet) sampleHetCount[i]++;
						else {
							if (hasRef) sampleHomoRefCount[i]++;
							if (hasAlt) sampleHomoAltCount[i]++;
						}

						if (hasAlt) sampleNonRefCount[i]++;
						if (hasRare) sampleRareAlleleCount[i]++;

						// Ts/Tv per sample (for biallelic SNPs)
						if (isBiallelicSnp && hasAlt) {
							if (isTransition) sampleTsCount[i]++;
							else              sampleTvCount[i]++;
						}
					}

					// Tajima's D
					if (isBiallelicSnp && siteAlleleCounts[0] > 0 && siteAlleleCounts[1] > 0) {
						numSegSites++;
						int n = totalGenotypedAlleles;
						if (n > 1) {
							double pi_i = (double)(siteAlleleCounts[0] * siteAlleleCounts[1]) / (n * (n - 1.0) / 2.0);
							piPerSiteSum += pi_i;
							totalAllelesForTajima += n;
							tajimaSiteCount++;
						}
					}
				} else {
					for (int i = 0; i < numSamples; i++) {
						if (parsedGenotypes[i] == null) sampleMissingCount[i]++;
					}
				}

				// HWE chi-square & Fisher Exact Test (biallelic diploid SNPs)
				if (isBiallelicSnp) {
					int nDip = n00 + n01 + n11;
					if (nDip >= 5) { // NGSEP uses > 3
						double pHwe = (2.0 * n00 + n01) / (2.0 * nDip);
						double qHwe = 1.0 - pHwe;

						// Clamping expected homozyogtes to avoid div-by-zero (like NGSEP)
						double e00 = Math.max(1.0, Math.min(nDip - 2.0, pHwe * pHwe * nDip));
						double e11 = Math.max(1.0, Math.min(nDip - 1.0 - e00, qHwe * qHwe * nDip));
						double e01 = Math.max(1.0, nDip - e00 - e11);

						double chi2 = (n00 - e00)*(n00 - e00)/e00
							+ (n01 - e01)*(n01 - e01)/e01
							+ (n11 - e11)*(n11 - e11)/e11;

						numHweTested++;
						chiSqSum += chi2;

						// Bins for histogram
						if      (chi2 < 2.0)    hweChiSqHistogram[0]++;
						else if (chi2 < 3.84)   hweChiSqHistogram[1]++;
						else if (chi2 < 7.0)    hweChiSqHistogram[2]++;
						else if (chi2 < 10.0)   hweChiSqHistogram[3]++;
						else                    hweChiSqHistogram[4]++;

						// Calculate Fisher exact test for P-value (like Wigginton 2005 / NGSEP)
						double fisherPValue = calculateHweFisherExactTest(n00, n01, n11);
						if (fisherPValue < 0.05) numHweViolations++;
					}
				}

				// Fst
				if (popAlleles != null && isBiallelicSnp) {
					int totRef = 0, totAlt = 0;
					for (int pi = 0; pi < numPops; pi++) {
						totRef += popAlleles[pi][0];
						totAlt += popAlleles[pi][1];
					}
					int totAll = totRef + totAlt;
					if (totAll > 0 && totRef > 0 && totAlt > 0) {
						double ptot = (double) totRef / totAll;
						double Ht   = 2.0 * ptot * (1.0 - ptot);
						for (int pi = 0; pi < numPops; pi++) {
							for (int pj = pi + 1; pj < numPops; pj++) {
								int nI = popAlleles[pi][0] + popAlleles[pi][1];
								int nJ = popAlleles[pj][0] + popAlleles[pj][1];
								if (nI == 0 || nJ == 0) continue;
								double pI = (double) popAlleles[pi][0] / nI;
								double pJ = (double) popAlleles[pj][0] / nJ;
								double Hs = (2.0*pI*(1-pI) + 2.0*pJ*(1-pJ)) / 2.0;
								double fst = (Ht > 0) ? (Ht - Hs) / Ht : 0;
								fstAccum[pi][pj][0] += Ht;
								fstAccum[pi][pj][1] += Hs;
								fstAccum[pi][pj][2]++;
							}
						}
					}
				}
			}
		}

		// ── Post-loop derived calculations ────────────────────────────────────
		meanEH     = ehCount > 0 ? ehSum / ehCount : 0.0;
		meanChiSq  = numHweTested > 0 ? chiSqSum / numHweTested : 0.0;

		sampleObsHet = new double[sampleNames.length];
		sampleFstat  = new double[sampleNames.length];
		for (int i = 0; i < sampleNames.length; i++) {
			double oh = sampleGenotypedCount[i] > 0
				? (double) sampleHetCount[i] / sampleGenotypedCount[i] : 0.0;
			sampleObsHet[i] = oh;
			sampleFstat[i]  = meanEH > 0 ? 1.0 - oh / meanEH : 0.0;
		}

		// ── Tajima's D ────────────────────────────────────────────────────────
		if (tajimaSiteCount > 0 && numSegSites > 1) {
			piHat = piPerSiteSum;
			int nEff = (int) Math.round((double) totalAllelesForTajima / tajimaSiteCount);
			if (nEff >= 2) {
				int S = numSegSites;
				double a1 = 0, a2 = 0;
				for (int i = 1; i < nEff; i++) {
					a1 += 1.0 / i;
					a2 += 1.0 / ((double)i * i);
				}
				thetaW = S / a1;
				double b1 = (nEff + 1.0) / (3.0 * (nEff - 1));
				double b2 = 2.0 * (nEff * nEff + nEff + 3.0) / (9.0 * nEff * (nEff - 1));
				double c1 = b1 - 1.0 / a1;
				double c2 = b2 - (nEff + 2.0) / (a1 * nEff) + a2 / (a1 * a1);
				double e1 = c1 / a1;
				double e2 = c2 / (a1 * a1 + a2);
				double varD = e1 * S + e2 * S * (S - 1);
				if (varD > 0) {
					tajimaD = (piHat - thetaW) / Math.sqrt(varD);
				}
			}
		}

		// ── Pairwise Fst ─────────────────────────────────────────────────────
		if (fstAccum != null) {
			pairwiseFst = new double[numPops][numPops];
			for (int pi = 0; pi < numPops; pi++) {
				for (int pj = pi + 1; pj < numPops; pj++) {
					double count = fstAccum[pi][pj][2];
					if (count > 0) {
						double meanHt = fstAccum[pi][pj][0] / count;
						double meanHs = fstAccum[pi][pj][1] / count;
						pairwiseFst[pi][pj] = meanHt > 0 ? (meanHt - meanHs) / meanHt : 0;
						pairwiseFst[pj][pi] = pairwiseFst[pi][pj];
					}
				}
			}
		}
	}

	// ── Fisher Exact Test for HWE (Wigginton et al. 2005) ────────────────────
	private double calculateHweFisherExactTest(int obsHom1, int obsHet, int obsHom2) {
		int n = obsHom1 + obsHet + obsHom2;
		int nA = 2 * obsHom1 + obsHet;
		int nB = 2 * obsHom2 + obsHet;

		int minorAlleleCount = Math.min(nA, nB);
		int majorAlleleCount = Math.max(nA, nB);

		int minHet = minorAlleleCount % 2;
		int maxHet = minorAlleleCount;

		double[] probs = new double[maxHet + 1];

		double[] logFacs = new double[n * 2 + 1];
		for (int i = 1; i <= n * 2; i++) logFacs[i] = logFacs[i - 1] + Math.log(i);

		double maxLogProb = -Double.MAX_VALUE;
		for (int h = minHet; h <= maxHet; h += 2) {
			int h1 = (minorAlleleCount - h) / 2;
			int h2 = n - h - h1;
			double logP = logFacs[n] + logFacs[minorAlleleCount] + logFacs[majorAlleleCount]
				- logFacs[2 * n] - logFacs[h1] - logFacs[h] - logFacs[h2] + h * Math.log(2.0);
			probs[h] = logP;
			if (logP > maxLogProb) maxLogProb = logP;
		}

		double sum = 0.0;
		for (int h = minHet; h <= maxHet; h += 2) {
			probs[h] = Math.exp(probs[h] - maxLogProb);
			sum += probs[h];
		}

		double pValue = 0.0;
		double obsProb = probs[obsHet];
		for (int h = minHet; h <= maxHet; h += 2) {
			if (probs[h] <= obsProb + 1e-9) {
				pValue += probs[h] / sum;
			}
		}
		return Math.min(1.0, pValue);
	}

	// ── Helpers ───────────────────────────────────────────────────────────────
	private boolean isTransition(String ref, String alt) {
		ref = ref.toUpperCase(); alt = alt.toUpperCase();
		return (ref.equals("A") && alt.equals("G")) ||
			   (ref.equals("G") && alt.equals("A")) ||
			   (ref.equals("C") && alt.equals("T")) ||
			   (ref.equals("T") && alt.equals("C"));
	}
}
