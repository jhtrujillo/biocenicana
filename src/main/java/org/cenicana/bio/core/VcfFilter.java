package org.cenicana.bio.core;
import org.cenicana.bio.utils.FileUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import org.cenicana.bio.io.VcfFastReader;

/**
 * A fast, streaming VCF filter that mimics NGSEP's VCFFilter.
 * Can filter by MAF, missingness, HWE p-value, and biallelic SNPs only.
 */
public class VcfFilter {

	private double minMaf = 0.0;
	private double maxMissingness = 1.0;
	private double minHwePValue = 0.0;
	private boolean onlyBiallelicSnps = false;
	private int ploidy = 2;
	private double minEh = 0.0;
	private int topN = -1;

	// Helper class for PriorityQueue to keep Top N polymorphic markers
	private static class VcfLineScore implements Comparable<VcfLineScore> {
		long order;
		String line;
		double eh;
		public VcfLineScore(long order, String line, double eh) {
			this.order = order;
			this.line = line;
			this.eh = eh;
		}
		@Override
		public int compareTo(VcfLineScore o) {
			return Double.compare(this.eh, o.eh); // Min-heap (lowest EH at the head)
		}
	}

	public void setMinMaf(double minMaf) { this.minMaf = minMaf; }
	public void setMaxMissingness(double maxMissingness) { this.maxMissingness = maxMissingness; }
	public void setMinHwePValue(double minHwePValue) { this.minHwePValue = minHwePValue; }
	public void setOnlyBiallelicSnps(boolean onlyBiallelicSnps) { this.onlyBiallelicSnps = onlyBiallelicSnps; }
	public void setPloidy(int ploidy) { this.ploidy = ploidy; }
	public void setMinEh(double minEh) { this.minEh = minEh; }
	public void setTopN(int topN) { this.topN = topN; }

	public void filter(String inputVcf, String outputVcf) throws IOException {
		String[] sampleNames = VcfFastReader.getSampleIds(inputVcf);
		int numSamples = sampleNames.length;

		int totalVariants = 0;
		int keptVariants = 0;

		PriorityQueue<VcfLineScore> topQueue = new PriorityQueue<>();

		try (BufferedReader br = new BufferedReader(new FileReader(inputVcf));
			 PrintWriter pw = new PrintWriter(new FileWriter(outputVcf))) {

			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					pw.println(line);
					if (line.startsWith("#CHROM")) {
						// Add a filter header before the CHROM line
						pw.println("##BioCenicana_VCFFilter=\"minMaf=" + minMaf + ",maxMissingness=" + maxMissingness + ",minHwePValue=" + minHwePValue + ",onlyBiallelicSnps=" + onlyBiallelicSnps + ",ploidy=" + ploidy + ",minEh=" + minEh + ",topN=" + topN + "\"");
					}
					continue;
				}

				String[] cols = line.split("\t");
				if (cols.length < 9 + numSamples) {
					pw.println(line); // Pass through malformed/short lines just in case
					continue;
				}

				totalVariants++;

				String ref = cols[3];
				String alt = cols[4];
				String[] altAllelesList = alt.equals(".") ? new String[0] : alt.split(",");
				int numAllelesAtSite = 1 + altAllelesList.length;

				boolean isBiallelicSnp = ref.length() == 1
					&& altAllelesList.length == 1
					&& altAllelesList[0].length() == 1
					&& !alt.equals(".");

				if (onlyBiallelicSnps && !isBiallelicSnp) {
					continue;
				}

				int siteMissing = 0;
				int totalGenotypedAlleles = 0;
				int[] siteAlleleCounts = new int[numAllelesAtSite];
				int n00 = 0, n01 = 0, n11 = 0; // HWE counts for diploid
				int[] polyDosageCounts = (ploidy > 2 && isBiallelicSnp) ? new int[ploidy + 1] : null;
				int nGenotypedPoly = 0;

				int formatAdIdx = -1;
				String[] formatFields = cols[8].split(":");
				for (int f = 0; f < formatFields.length; f++) {
					if (formatFields[f].equals("AD") || formatFields[f].equals("BSDP")) {
						formatAdIdx = f; break;
					}
				}

				for (int i = 0; i < numSamples; i++) {
					String[] gData = cols[9 + i].split(":");
					String gt = gData[0];

					if (ploidy > 2 && formatAdIdx != -1 && isBiallelicSnp) {
						// Polyploid mode: rely on AD/BSDP to estimate dosage
						boolean dosageEstimated = false;
						if (gData.length > formatAdIdx && !gData[formatAdIdx].equals(".")) {
							String[] adStrs = gData[formatAdIdx].split(",");
							if (adStrs.length >= 2) {
								try {
									double refReads = Double.parseDouble(adStrs[0]);
									double altReads = Double.parseDouble(adStrs[1]);
									double totalReads = refReads + altReads;
									if (totalReads >= 5) { // Minimum reads to confidently estimate dosage
										int dosage = (int) Math.round((altReads / totalReads) * ploidy);
										if (dosage >= 0 && dosage <= ploidy) {
											polyDosageCounts[dosage]++;
											siteAlleleCounts[0] += (ploidy - dosage);
											siteAlleleCounts[1] += dosage;
											totalGenotypedAlleles += ploidy;
											nGenotypedPoly++;
											dosageEstimated = true;
										}
									}
								} catch (NumberFormatException ignored) {}
							}
						}
						if (!dosageEstimated) {
							siteMissing++;
						}
					} else {
						// Standard diploid/explicit genotype mode
						if (gt.equals(".") || gt.startsWith("./.") || gt.startsWith("./")) {
							siteMissing++;
						} else {
							String[] alleles = gt.split("[/|]");
							boolean valid = true;
							int[] numAlleles = new int[alleles.length];
							for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
								if (alleles[aIdx].equals(".")) {
									valid = false;
									break;
								}
								try {
									numAlleles[aIdx] = Integer.parseInt(alleles[aIdx]);
									if (numAlleles[aIdx] >= numAllelesAtSite) valid = false;
								} catch (NumberFormatException e) {
									valid = false;
								}
							}

							if (valid) {
								for (int a : numAlleles) {
									siteAlleleCounts[a]++;
									totalGenotypedAlleles++;
								}
								if (isBiallelicSnp && numAlleles.length == 2) {
									if (numAlleles[0] == 0 && numAlleles[1] == 0) n00++;
									else if (numAlleles[0] == 1 && numAlleles[1] == 1) n11++;
									else n01++;
								}
							} else {
								siteMissing++;
							}
						}
					}
				}

				// --- Filtering Checks ---

				// 1. Missingness filter
				double missingFreq = (double) siteMissing / numSamples;
				if (missingFreq > maxMissingness) {
					continue;
				}

				if (totalGenotypedAlleles > 0) {
					// 2. MAF filter
					int minAc = Integer.MAX_VALUE;
					for (int c : siteAlleleCounts) {
						if (c > 0 && c < minAc) {
							minAc = c;
						}
					}
					double maf = (double) minAc / totalGenotypedAlleles;
					if (maf < minMaf) {
						continue;
					}

					// 3. HWE Filter
					if (isBiallelicSnp && minHwePValue > 0.0) {
						if (ploidy > 2 && polyDosageCounts != null && nGenotypedPoly >= 5) {
							// Polyploid Polysomic HWE Chi-square test
							double totalAltDosage = 0;
							for (int d = 0; d <= ploidy; d++) totalAltDosage += polyDosageCounts[d] * d;
							double q = totalAltDosage / (nGenotypedPoly * ploidy);
							double p = 1.0 - q;

							double chi2 = 0.0;
							double[] logFacs = new double[ploidy + 1];
							for(int j=1; j<=ploidy; j++) logFacs[j] = logFacs[j-1] + Math.log(j);

							for (int d = 0; d <= ploidy; d++) {
								double logComb = logFacs[ploidy] - logFacs[d] - logFacs[ploidy - d];
								double expectedProb = Math.exp(logComb + (ploidy - d)*Math.log(Math.max(1e-10, p)) + d*Math.log(Math.max(1e-10, q)));
								double expectedCount = Math.max(1e-5, nGenotypedPoly * expectedProb);
								double observedCount = polyDosageCounts[d];
								chi2 += (observedCount - expectedCount) * (observedCount - expectedCount) / expectedCount;
							}

							double pValue = chiSquarePValue(chi2, ploidy - 1);
							if (pValue < minHwePValue) {
								continue;
							}
						} else if (ploidy <= 2) {
							// Diploid exact test
							int nDip = n00 + n01 + n11;
							if (nDip >= 5) {
								double fisherPValue = calculateHweFisherExactTest(n00, n01, n11);
								if (fisherPValue < minHwePValue) {
									continue;
								}
							}
						}
					}
					// 4. Expected Heterozygosity (EH) filter
					double eh = 0.0;
					for (int c : siteAlleleCounts) {
						if (c > 0) {
							double freq = (double) c / totalGenotypedAlleles;
							eh += freq * freq;
						}
					}
					eh = 1.0 - eh;

					if (eh < minEh) {
						continue;
					}

					// If we are keeping Top N, add to PriorityQueue instead of writing immediately
					if (topN > 0) {
						if (topQueue.size() < topN) {
							topQueue.add(new VcfLineScore(totalVariants, line, eh));
						} else if (eh > topQueue.peek().eh) {
							topQueue.poll();
							topQueue.add(new VcfLineScore(totalVariants, line, eh));
						}
						continue; // Skip the direct writing below
					}

				} else {
					continue; // Filter out sites with 100% missing data
				}

				// If it passed all filters and we are not in Top N mode, keep it
				pw.println(line);
				keptVariants++;
			}

			// If Top N mode is active, write the buffered markers in their original genomic order
			if (topN > 0) {
				List<VcfLineScore> sortedList = new ArrayList<>(topQueue);
				sortedList.sort(Comparator.comparingLong(s -> s.order));
				for (VcfLineScore s : sortedList) {
					pw.println(s.line);
					keptVariants++;
				}
			}
		}

		System.out.println("[VcfFilter] Processing complete.");
		System.out.println("[VcfFilter] Total variants processed: " + totalVariants);
		System.out.println("[VcfFilter] Variants kept: " + keptVariants);
		System.out.println("[VcfFilter] Variants removed: " + (totalVariants - keptVariants));
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

	private double chiSquarePValue(double x, int df) {
		if (x <= 0 || df < 1) return 1.0;
		double dfD = df;
		// Wilson-Hilferty transformation for Chi-Square to Normal distribution
		double Z = (Math.pow(x / dfD, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * dfD))) / Math.sqrt(2.0 / (9.0 * dfD));
		return 0.5 * erfc(Z / Math.sqrt(2.0));
	}

	private double erfc(double x) {
		double z = Math.abs(x);
		double t = 1.0 / (1.0 + 0.5 * z);
		double ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
		return x >= 0 ? ans : 2.0 - ans;
	}
}
