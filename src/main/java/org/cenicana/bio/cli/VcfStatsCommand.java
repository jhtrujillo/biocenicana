package org.cenicana.bio.cli;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.concurrent.Callable;
import org.cenicana.bio.VcfStatisticsCalculator;
import org.cenicana.bio.io.HtmlDashboardGenerator;

@Command(name = "vcf-stats", description = "Generates a high-quality HTML statistics dashboard from a VCF file (supports NGSEP, GATK, FreeBayes)", mixinStandardHelpOptions = true)
public class VcfStatsCommand implements Callable<Integer> {

	@Option(names = { "-v", "--vcf" }, required = true, description = "Path to the input VCF file")
	private String vcfFile;

	@Option(names = { "-o", "--output" }, required = true, description = "Path to the output HTML dashboard file (e.g. report.html)")
	private String outputFile;

	@Override
	public Integer call() throws Exception {
		System.out.println("Analyzing VCF file: " + vcfFile);
		
		VcfStatisticsCalculator calc = new VcfStatisticsCalculator();
		calc.calculate(vcfFile);
		
		System.out.println("Analysis complete. Found " + (calc.numSnps + calc.numIndels) + " total variants.");
		System.out.println("Generating High-Quality HTML Dashboard at: " + outputFile);

		HtmlDashboardGenerator.generateReport(calc, outputFile);

		System.out.println("Success! Open " + outputFile + " in any web browser.");
		return 0;
	}
}
