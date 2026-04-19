package org.cenicana.bio.cli;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.io.File;
import java.util.concurrent.Callable;
import org.cenicana.bio.VcfStatisticsCalculator;
import org.cenicana.bio.io.HtmlDashboardGenerator;
import org.cenicana.bio.io.TsvStatsExporter;

@Command(name = "vcf-stats",
	description = "Generates a HTML dashboard and TSV summary from a VCF file. "
		+ "Results are saved inside a folder named after the value given to -o. "
		+ "(supports NGSEP, GATK, FreeBayes)",
	mixinStandardHelpOptions = true)
public class VcfStatsCommand implements Callable<Integer> {

	@Option(names = { "-v", "--vcf" }, required = true, description = "Path to the input VCF file")
	private String vcfFile;

	@Option(names = { "-o", "--output" }, required = true,
		description = "Base name for the output folder and files. "
			+ "E.g. '-o test' creates test/test.html and test/test.tsv")
	private String baseName;

	@Override
	public Integer call() throws Exception {
		System.out.println("[biocenicana] Analyzing VCF: " + vcfFile);

		// Create output directory
		File outDir = new File(baseName);
		if (!outDir.exists()) {
			outDir.mkdirs();
			System.out.println("[biocenicana] Created output folder: " + outDir.getAbsolutePath());
		}

		// Run statistics
		VcfStatisticsCalculator calc = new VcfStatisticsCalculator();
		calc.calculate(vcfFile);
		System.out.println("[biocenicana] Analysis complete. Found "
			+ (calc.numSnps + calc.numIndels) + " variants across "
			+ calc.sampleNames.length + " samples.");

		// Build file paths: <baseName>/<baseName>.html and <baseName>/<baseName>.tsv
		String htmlPath = outDir.getPath() + File.separator + baseName + ".html";
		String tsvPath  = outDir.getPath() + File.separator + baseName + ".tsv";

		// Export HTML dashboard
		HtmlDashboardGenerator.generateReport(calc, htmlPath);
		System.out.println("[biocenicana] HTML Dashboard -> " + htmlPath);

		// Export TSV summary
		TsvStatsExporter.exportSampleStats(calc, tsvPath);
		System.out.println("[biocenicana] TSV Summary    -> " + tsvPath);

		System.out.println("\n✅ Done! Open the following file in your browser:");
		System.out.println("   " + htmlPath);
		return 0;
	}
}
