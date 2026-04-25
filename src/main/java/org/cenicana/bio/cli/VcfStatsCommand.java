package org.cenicana.bio.cli;
import org.cenicana.bio.utils.FileUtils;
import org.cenicana.bio.core.JoinMapCpFormat;
import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.VcfStatisticsCalculator;
import org.cenicana.bio.core.VcfFilter;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.io.File;
import java.util.concurrent.Callable;
import org.cenicana.bio.io.HtmlDashboardGenerator;
import org.cenicana.bio.io.TsvStatsExporter;

@Command(name = "vcf-stats",
	description = "Generates a HTML dashboard + TSV summary from a VCF file. "
		+ "Supports NGSEP, GATK and FreeBayes. "
		+ "Outputs are saved inside a folder named after the -o value.",
	mixinStandardHelpOptions = true)
public class VcfStatsCommand implements Callable<Integer> {

	@Option(names = { "-v", "--vcf" }, required = true,
		description = "Path to the input VCF file.")
	private String vcfFile;

	@Option(names = { "-o", "--output" }, required = true,
		description = "Base name for the output folder. "
			+ "E.g. '-o test' creates test/test.html and test/test.tsv.")
	private String baseName;

	@Option(names = { "-p", "--ploidy" }, description = "Ploidy level of the population (optional for stats).")
	private Integer ploidy;

	@Option(names = { "--popmap" }, required = false,
		description = "Optional: tab-delimited file with sample populations. "
			+ "Two columns: SAMPLE_ID and POPULATION_ID. "
			+ "When provided, pairwise Fst is calculated between populations.")
	private String popFile;

	@Option(names = { "-t", "--threads" }, description = "Number of threads for parallel processing.", defaultValue = "-1")
	private int threads;

	@Override
	public Integer call() throws Exception {
		System.out.println("[biocenicana] Analyzing VCF: " + vcfFile);

		// Create output directory
		File outDir = new File(baseName);
		if (!outDir.exists()) {
			outDir.mkdirs();
			System.out.println("[biocenicana] Created output folder: " + outDir.getAbsolutePath());
		}

		// Initialize calculator
		VcfStatisticsCalculator calc = new VcfStatisticsCalculator();

		// Load population file if provided (enables Fst)
		if (popFile != null) {
			System.out.println("[biocenicana] Loading population assignments from: " + popFile);
			calc.loadPopulationMap(popFile);
			System.out.println("[biocenicana] Populations found: "
				+ calc.populationNames.length + " → " + String.join(", ", calc.populationNames));
		}

		int numThreads = (threads <= 0) ? Runtime.getRuntime().availableProcessors() : threads;
		calc.calculate(vcfFile, numThreads);
		System.out.println("[biocenicana] Done. Variants: "
			+ (calc.numSnps + calc.numIndels)
			+ "  |  Samples: " + calc.sampleNames.length
			+ "  |  Seg. sites: " + calc.numSegSites);

		if (!Double.isNaN(calc.tajimaD)) {
			System.out.printf("[biocenicana] Tajima's D = %.4f%n", calc.tajimaD);
		}
		if (calc.numHweTested > 0) {
			double pctHwe = 100.0 * calc.numHweViolations / calc.numHweTested;
			System.out.printf("[biocenicana] HWE violations: %d / %d sites (%.1f%%)%n",
				calc.numHweViolations, calc.numHweTested, pctHwe);
		}

		// Build output file paths
		String simpleName = new File(baseName).getName();
		String htmlPath = outDir.getPath() + File.separator + simpleName + ".html";
		String tsvPath  = outDir.getPath() + File.separator + simpleName + ".tsv";

		HtmlDashboardGenerator.generateReport(calc, htmlPath);
		System.out.println("[biocenicana] HTML Dashboard → " + htmlPath);

		TsvStatsExporter.exportSampleStats(calc, tsvPath);
		System.out.println("[biocenicana] TSV Summary    → " + tsvPath);

		System.out.println("\n✅  Open in your browser: " + htmlPath);
		return 0;
	}
}
