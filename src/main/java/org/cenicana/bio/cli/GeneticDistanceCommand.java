package org.cenicana.bio.cli;

import java.util.concurrent.Callable;

import org.cenicana.bio.AlleleDosageCalculator;

import picocli.CommandLine.Command;
import picocli.CommandLine.ITypeConverter;
import picocli.CommandLine.Option;
import picocli.CommandLine.TypeConversionException;

@Command(
	name = "genetic-distance",
	description = "Calculates the N x N genetic distance matrix between individuals from a VCF file.",
	mixinStandardHelpOptions = true,
	version = "genetic-distance 1.0"
)
public class GeneticDistanceCommand implements Callable<Integer> {

	// ── Enums & Converters ─────────────────────────────────────────────────────

	enum CallerType {
		NGSEP("ngsep"),
		GATK("gatk"),
		FREEBAYES("freebayes"),
		AUTO("auto");

		private final String stringValue;

		CallerType(String stringValue) {
			this.stringValue = stringValue;
		}

		@Override
		public String toString() {
			return stringValue;
		}
	}

	static class CallerConverter implements ITypeConverter<CallerType> {
		@Override
		public CallerType convert(String value) throws TypeConversionException {
			for (CallerType c : CallerType.values()) {
				if (c.stringValue.equalsIgnoreCase(value)) {
					return c;
				}
			}
			throw new TypeConversionException(
				"Invalid value: " + value + ". Allowed values: ngsep, gatk, freebayes, auto"
			);
		}
	}

	// ── Options ────────────────────────────────────────────────────────────────

	@Option(names = {"-v", "--vcf"},
		description = "Path to the VCF file.",
		required = true)
	private String vcfFile;

	@Option(names = {"-p", "--ploidy"},
		description = "Ploidy level of the organism (e.g., 10 for sugarcane).",
		required = true)
	private int ploidy;

	@Option(names = {"-c", "--caller"},
		description = "Variant caller used (ngsep, gatk, freebayes, auto). Default: ${DEFAULT-VALUE}.",
		defaultValue = "auto",
		converter = CallerConverter.class)
	private CallerType caller;

	@Option(names = {"-md", "--min-depth"},
		description = "Minimum total read depth (Ref + Alt) required to trust a genotype call. Lower depth calls will be treated as missing data. Default: ${DEFAULT-VALUE}.",
		defaultValue = "0")
	private int minDepth;

	@Option(names = {"-ar", "--adaptive-rounding"},
		description = "Enable Machine Learning (K-Means 1D) adaptive rounding to overcome reference bias before calculating distances.",
		defaultValue = "false")
	private boolean adaptiveRounding;

	// ── Execution ──────────────────────────────────────────────────────────────

	@Override
	public Integer call() throws Exception {
		try {
			AlleleDosageCalculator calculator = new AlleleDosageCalculator();
			calculator.computeAndPrintDistanceMatrix(vcfFile, ploidy, caller.toString(), minDepth, adaptiveRounding);
			return 0;

		} catch (Exception e) {
			System.err.println("Error running genetic-distance: " + e.getMessage());
			e.printStackTrace();
			return 1;
		}
	}
}
