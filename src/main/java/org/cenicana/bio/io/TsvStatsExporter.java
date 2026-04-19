package org.cenicana.bio.io;

import org.cenicana.bio.VcfStatisticsCalculator;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TsvStatsExporter {

	public static void exportSampleStats(VcfStatisticsCalculator stats, String outputPath) throws IOException {
		int numVariants = stats.numSnps + stats.numIndels;
		if (numVariants == 0) numVariants = 1;

		try (PrintWriter writer = new PrintWriter(new FileWriter(outputPath))) {
			// Header
			writer.println("## BioCenicana VCF Statistics - Per Sample Summary");
			writer.println("## Total Variants: " + (stats.numSnps + stats.numIndels));
			writer.println("## Total SNPs: " + stats.numSnps);
			writer.println("## Total InDels: " + stats.numIndels);
			double tsTv = stats.numTransversions == 0 ? 0 : (double) stats.numTransitions / stats.numTransversions;
			writer.printf("## Ts/Tv Ratio: %.3f%n", tsTv);
			writer.println("#");
			writer.println("SAMPLE_ID\tAVG_DEPTH\tMISSING_COUNT\tMISSING_RATE_PCT");

			for (int i = 0; i < stats.sampleNames.length; i++) {
				double avgDepth = stats.sampleDepthCount[i] == 0 ? 0.0
						: (double) stats.sampleTotalDepth[i] / stats.sampleDepthCount[i];
				double missingPct = ((double) stats.sampleMissingCount[i] / numVariants) * 100.0;

				writer.printf("%s\t%.2f\t%d\t%.2f%n",
						stats.sampleNames[i],
						avgDepth,
						stats.sampleMissingCount[i],
						missingPct);
			}
		}
	}
}
