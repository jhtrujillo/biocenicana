package org.cenicana.bio.cli;
import org.cenicana.bio.utils.FileUtils;
import org.cenicana.bio.core.JoinMapCpFormat;
import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.VcfStatisticsCalculator;
import org.cenicana.bio.core.VcfFilter;

import java.io.File;
import java.util.concurrent.Callable;


import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command(name = "vcf-filter", description = "Filter VCF file by MAF, missingness, and HWE p-value (similar to NGSEP VCFFilter)")
public class VcfFilterCommand implements Callable<Integer> {

	@Option(names = {"-v", "--vcf"}, required = true, description = "Input VCF file")
	private String vcfFile;

	@Option(names = {"-o", "--output"}, required = true, description = "Output VCF file")
	private String outputFile;

	@Option(names = {"--minMAF"}, description = "Minimum Minor Allele Frequency (MAF)", defaultValue = "0.0")
	private double minMaf;

	@Option(names = {"-m", "--min-samples"}, description = "Minimum number of samples genotyped to keep the variant (similar to NGSEP -m)", defaultValue = "0")
	private int minSamples;

	@Option(names = {"-x", "--max-missingness"}, description = "Maximum fraction of missing genotypes (0.0 to 1.0). If -m is provided, this is ignored.", defaultValue = "1.0")
	private double maxMissingness;

	@Option(names = {"-e", "--min-hwe-pvalue"}, description = "Minimum HWE p-value (Fisher exact test). Use 0.05 to remove SNPs with HWE violation", defaultValue = "0.0")
	private double minHwePValue;

	@Option(names = {"-b", "--biallelic-only"}, description = "Keep only biallelic SNPs (remove indels and multiallelics)", defaultValue = "false")
	private boolean onlyBiallelicSnps;

	@Option(names = {"-p", "--ploidy"}, description = "Ploidy level of the organism. If > 2, HWE and MAF will be calculated using Allele Dosages (AD/BSDP) assuming polysomic inheritance.", defaultValue = "2")
	private int ploidy;

	@Option(names = {"--min-eh"}, description = "Minimum Expected Heterozygosity (EH) score. Ranges from 0.0 to 1.0.", defaultValue = "0.0")
	private double minEh;

	@Option(names = {"--top-n"}, description = "Keep only the Top N most polymorphic markers based on Expected Heterozygosity (EH).", defaultValue = "-1")
	private int topN;

	@Option(names = {"-d", "--min-rd"}, description = "Minimum read depth (Ref + Alt) required to trust a genotype call in polyploid mode. Default: ${DEFAULT-VALUE}.", defaultValue = "5")
	private int minDepth;

	@Option(names = {"-s", "--strict"}, description = "Strict mode: discard genotypes marked as ./. even if they have enough read counts.", defaultValue = "false")
	private boolean strictGenotypes;

	@Option(names = {"-t", "--threads"}, description = "Number of threads to use for parallel processing. Default: number of available processors.", defaultValue = "-1")
	private int threads;

	@Override
	public Integer call() throws Exception {
		File f = new File(vcfFile);
		if (!f.exists()) {
			System.err.println("Error: VCF file not found: " + vcfFile);
			return 1;
		}

		System.out.println("=================================================");
		System.out.println("BioCenicana VCF Filter");
		System.out.println("=================================================");
		System.out.println("Input:  " + vcfFile);
		System.out.println("Output: " + outputFile);
		System.out.println("Filters applied:");
		if (minMaf > 0) System.out.println(" - Minimum MAF (--minMAF): " + minMaf);
		if (minSamples > 0) System.out.println(" - Min Genotyped Samples (-m): " + minSamples);
		else if (maxMissingness < 1.0) System.out.println(" - Max Missingness (-x): " + maxMissingness);
		if (minHwePValue > 0) System.out.println(" - Min HWE p-value: " + minHwePValue);
		if (onlyBiallelicSnps) System.out.println(" - Biallelic SNPs only: true");
		System.out.println(" - Ploidy assumption: " + ploidy);
		if (minEh > 0) System.out.println(" - Minimum EH: " + minEh);
		if (topN > 0) System.out.println(" - Keep Top N polymorphic: " + topN);
		if (ploidy > 2) System.out.println(" - Minimum Read Depth: " + minDepth);
		if (strictGenotypes) System.out.println(" - Strict Mode (discard ./.): true");
		System.out.println("=================================================\n");

		VcfFilter filter = new VcfFilter();
		filter.setMinMaf(minMaf);
		
		// Convert minSamples to missingness if provided
		if (minSamples > 0) {
			int totalSamples = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile).length;
			double missingRate = (double)(totalSamples - minSamples) / totalSamples;
			filter.setMaxMissingness(missingRate);
		} else {
			filter.setMaxMissingness(maxMissingness);
		}
		
		filter.setMinHwePValue(minHwePValue);
		filter.setOnlyBiallelicSnps(onlyBiallelicSnps);
		filter.setPloidy(ploidy);
		filter.setMinEh(minEh);
		filter.setTopN(topN);
		filter.setMinDepth(minDepth);
		filter.setStrictGenotypes(strictGenotypes);

		int numThreads = (threads <= 0) ? Runtime.getRuntime().availableProcessors() : threads;
		filter.filter(vcfFile, outputFile, numThreads);

		return 0;
	}
}
