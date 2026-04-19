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

		try (PrintWriter writer = new PrintWriter(new FileWriter(outputPath))) {
			// Global summary header
			writer.println("## BioCenicana VCF Statistics - Per Sample Summary");
			writer.printf("## Total Variants: %d  (SNPs: %d | InDels: %d)%n",
				numVariants, stats.numSnps, stats.numIndels);
			writer.printf("## Transitions: %d  |  Transversions: %d%n",
				stats.numTransitions, stats.numTransversions);
			double tsTv = stats.numTransversions == 0 ? 0
				: (double) stats.numTransitions / stats.numTransversions;
			writer.printf("## Ts/Tv Ratio: %.3f%n", tsTv);
			writer.printf("## Mean Expected Heterozygosity (EH): %.4f%n", stats.meanEH);
			writer.println("#");

			// Variant density per chromosome
			writer.println("## -- Variant Density per Chromosome --");
			for (Map.Entry<String, Integer> entry : stats.variantsPerChromosome.entrySet()) {
				writer.printf("##   %s: %d variants%n", entry.getKey(), entry.getValue());
			}
			writer.println("#");

			// Per-sample table
			writer.println("SAMPLE_ID\tAVG_DEPTH\tMISSING_COUNT\tMISSING_RATE_PCT\tGENOTYPED_COUNT\tHET_COUNT\tOBS_HET\tF_STATISTIC");

			for (int i = 0; i < stats.sampleNames.length; i++) {
				double avgDepth = stats.sampleDepthCount[i] == 0 ? 0.0
					: (double) stats.sampleTotalDepth[i] / stats.sampleDepthCount[i];
				double missingPct = ((double) stats.sampleMissingCount[i] / numVariants) * 100.0;

				writer.printf("%s\t%.2f\t%d\t%.2f\t%d\t%d\t%.4f\t%.4f%n",
					stats.sampleNames[i],
					avgDepth,
					stats.sampleMissingCount[i],
					missingPct,
					stats.sampleGenotypedCount[i],
					stats.sampleHetCount[i],
					stats.sampleObsHet[i],
					stats.sampleFstat[i]);
			}
		}
	}
}
