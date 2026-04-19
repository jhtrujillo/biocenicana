package org.cenicana.bio.cli;

import java.io.File;
import java.util.concurrent.Callable;

import org.cenicana.bio.VcfFilter;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command(name = "vcf-filter", description = "Filter VCF file by MAF, missingness, and HWE p-value (similar to NGSEP VCFFilter)")
public class VcfFilterCommand implements Callable<Integer> {

	@Option(names = {"-v", "--vcf"}, required = true, description = "Input VCF file")
	private String vcfFile;

	@Option(names = {"-o", "--output"}, required = true, description = "Output VCF file")
	private String outputFile;

	@Option(names = {"-m", "--min-maf"}, description = "Minimum Minor Allele Frequency (MAF)", defaultValue = "0.0")
	private double minMaf;

	@Option(names = {"-x", "--max-missingness"}, description = "Maximum fraction of missing genotypes (0.0 to 1.0)", defaultValue = "1.0")
	private double maxMissingness;

	@Option(names = {"-e", "--min-hwe-pvalue"}, description = "Minimum HWE p-value (Fisher exact test). Use 0.05 to remove SNPs with HWE violation", defaultValue = "0.0")
	private double minHwePValue;

	@Option(names = {"-b", "--biallelic-only"}, description = "Keep only biallelic SNPs (remove indels and multiallelics)", defaultValue = "false")
	private boolean onlyBiallelicSnps;

	@Option(names = {"-p", "--ploidy"}, description = "Ploidy level of the organism. If > 2, HWE and MAF will be calculated using Allele Dosages (AD/BSDP) assuming polysomic inheritance.", defaultValue = "2")
	private int ploidy;

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
		if (minMaf > 0) System.out.println(" - Minimum MAF: " + minMaf);
		if (maxMissingness < 1.0) System.out.println(" - Max Missingness: " + maxMissingness);
		if (minHwePValue > 0) System.out.println(" - Min HWE p-value: " + minHwePValue);
		if (onlyBiallelicSnps) System.out.println(" - Biallelic SNPs only: true");
		System.out.println(" - Ploidy assumption: " + ploidy);
		System.out.println("=================================================\n");

		VcfFilter filter = new VcfFilter();
		filter.setMinMaf(minMaf);
		filter.setMaxMissingness(maxMissingness);
		filter.setMinHwePValue(minHwePValue);
		filter.setOnlyBiallelicSnps(onlyBiallelicSnps);
		filter.setPloidy(ploidy);

		filter.filter(vcfFile, outputFile);

		return 0;
	}
}
