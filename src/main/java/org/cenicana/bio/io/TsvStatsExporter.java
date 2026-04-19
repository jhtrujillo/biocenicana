package org.cenicana.bio.io;

import org.cenicana.bio.VcfStatisticsCalculator;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

public class TsvStatsExporter {

	public static void exportSampleStats(VcfStatisticsCalculator stats, String outputPath) throws IOException {
		int numVariants = stats.numSnps + stats.numIndels;
		if (numVariants == 0) numVariants = 1;

		try (PrintWriter w = new PrintWriter(new FileWriter(outputPath))) {

			// ── SECTION 1: Global summary ─────────────────────────────────────
			w.println("## ========================================================");
			w.println("## BioCenicana VCF Statistics Report");
			w.println("## ========================================================");
			w.printf("## Total Variants: %d  (SNPs: %d | InDels: %d)%n",
				numVariants, stats.numSnps, stats.numIndels);
			w.printf("## Transitions: %d  |  Transversions: %d%n",
				stats.numTransitions, stats.numTransversions);
			double tsTv = stats.numTransversions == 0 ? 0
				: (double) stats.numTransitions / stats.numTransversions;
			w.printf("## Ts/Tv Ratio: %.3f%n", tsTv);
			w.printf("## Mean Expected Heterozygosity (EH): %.4f%n", stats.meanEH);
			w.println("#");

			// ── SECTION 2: Level 2 stats ──────────────────────────────────────
			w.println("## -- AN: Allele Number Distribution --");
			w.printf("##   Monomorphic sites: %d%n", stats.numMonomorphic);
			w.printf("##   Biallelic sites:   %d%n", stats.numBiallelic);
			w.printf("##   Multiallelic sites:%d%n", stats.numMultiallelic);
			w.println("#");

			w.println("## -- HWE Chi-square Test --");
			w.printf("##   Sites tested: %d%n", stats.numHweTested);
			w.printf("##   Sites violating HWE (chi2 > 3.84, p<0.05): %d (%.1f%%)%n",
				stats.numHweViolations,
				stats.numHweTested > 0 ? 100.0 * stats.numHweViolations / stats.numHweTested : 0.0);
			w.printf("##   Mean chi-square: %.3f%n", stats.meanChiSq);
			w.println("#");

			w.println("## -- Tajima's D --");
			w.printf("##   Segregating sites (S): %d%n", stats.numSegSites);
			w.printf("##   pi (avg pairwise diffs): %.5f%n", stats.piHat);
			w.printf("##   Watterson theta_W: %.5f%n", stats.thetaW);
			w.printf("##   Tajima's D: %s%n",
				Double.isNaN(stats.tajimaD) ? "N/A (insufficient data)"
				: String.format("%.4f", stats.tajimaD));
			if (!Double.isNaN(stats.tajimaD)) {
				String interpretation;
				if      (stats.tajimaD < -2.0) interpretation = "Strong negative selection or population expansion";
				else if (stats.tajimaD < 0)    interpretation = "Slight tendency toward negative selection / expansion";
				else if (stats.tajimaD < 2.0)  interpretation = "Neutral evolution (consistent with drift)";
				else                            interpretation = "Balancing selection or population contraction";
				w.printf("##   Interpretation: %s%n", interpretation);
			}
			w.println("#");

			// ── SECTION 3: Fst (if populations were provided) ────────────────
			if (stats.pairwiseFst != null && stats.populationNames != null) {
				w.println("## -- Pairwise Fst --");
				int np = stats.populationNames.length;
				for (int i = 0; i < np; i++) {
					for (int j = i + 1; j < np; j++) {
						w.printf("##   Fst %s vs %s: %.4f%n",
							stats.populationNames[i],
							stats.populationNames[j],
							stats.pairwiseFst[i][j]);
					}
				}
				w.println("#");
			}

			// ── SECTION 4: Density per chromosome ────────────────────────────
			w.println("## -- Variant Density per Chromosome --");
			for (Map.Entry<String, Integer> entry : stats.variantsPerChromosome.entrySet()) {
				w.printf("##   %s: %d variants%n", entry.getKey(), entry.getValue());
			}
			w.println("#");

			// ── SECTION 5: Per-sample table ───────────────────────────────────
			w.println("SAMPLE_ID\tAVG_DEPTH\tMISSING_COUNT\tMISSING_RATE_PCT\tGENOTYPED_COUNT\tHOMO_REF_COUNT\tHOMO_ALT_COUNT\tHET_COUNT\tNON_REF_COUNT\tRARE_ALLELE_COUNT\tTS_COUNT\tTV_COUNT\tTS_TV_RATIO\tOBS_HET\tF_STATISTIC");

			for (int i = 0; i < stats.sampleNames.length; i++) {
				double avgDepth = stats.sampleDepthCount[i] == 0 ? 0.0
					: (double) stats.sampleTotalDepth[i] / stats.sampleDepthCount[i];
				double missingPct = ((double) stats.sampleMissingCount[i] / numVariants) * 100.0;
				double sampleTsTv = stats.sampleTvCount[i] == 0 ? 0.0 : (double) stats.sampleTsCount[i] / stats.sampleTvCount[i];

				w.printf("%s\t%.2f\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.4f\t%.4f%n",
					stats.sampleNames[i],
					avgDepth,
					stats.sampleMissingCount[i],
					missingPct,
					stats.sampleGenotypedCount[i],
					stats.sampleHomoRefCount[i],
					stats.sampleHomoAltCount[i],
					stats.sampleHetCount[i],
					stats.sampleNonRefCount[i],
					stats.sampleRareAlleleCount[i],
					stats.sampleTsCount[i],
					stats.sampleTvCount[i],
					sampleTsTv,
					stats.sampleObsHet[i],
					stats.sampleFstat[i]);
			}
		}
	}
}
