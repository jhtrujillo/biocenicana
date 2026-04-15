package huellamolecular.cli;

import huellamolecular.GeneDosis;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.ITypeConverter;
import picocli.CommandLine.TypeConversionException;
import java.util.concurrent.Callable;

@Command(
	name = "allele-dosage",
	mixinStandardHelpOptions = true,
	description = "Computes an allele dosage matrix from a VCF file for polyploid organisms.",
	sortOptions = false
)
public class AlleleDosageCommand implements Callable<Integer> {

	// ── Imputation method ──────────────────────────────────────────────────────

	public enum ImputationMethod {
		bsdp      ("bsdp",      "Uses allele counts from VCF. Missing genotypes left as -1."),
		mode      ("mode",      "Imputes missing values with the mode dosage of the SNP."),
		mean      ("mean",      "Imputes missing values with the mean dosage of the SNP."),
		bsdp_mode ("bsdp-mode", "Counts first; remaining missing values imputed with the mode."),
		bsdp_mean ("bsdp-mean", "Counts first; remaining missing values imputed with the mean.");

		private final String cliName;

		ImputationMethod(String cliName, String description) {
			this.cliName = cliName;
		}

		public String toInternalKey() { return cliName; }

		@Override
		public String toString() { return cliName; }
	}

	static class ImputationConverter implements ITypeConverter<ImputationMethod> {
		@Override
		public ImputationMethod convert(String value) throws TypeConversionException {
			String normalized = value.replace('-', '_').toLowerCase();
			for (ImputationMethod m : ImputationMethod.values()) {
				if (m.name().equals(normalized)) return m;
			}
			throw new TypeConversionException(
				"Invalid imputation method '" + value + "'. Valid: bsdp, mode, mean, bsdp-mode, bsdp-mean");
		}
	}

	// ── Variant caller type ────────────────────────────────────────────────────

	public enum CallerType {
		auto      ("auto",      "Auto-detect the allele depth field from FORMAT column."),
		ngsep     ("ngsep",     "Use NGSEP's BSDP field [countAlt, countRef, ...]."),
		gatk      ("gatk",      "Use GATK's AD field [countRef, countAlt]."),
		freebayes ("freebayes", "Use FreeBayes's RO (ref) and AO (alt) fields.");

		private final String cliName;

		CallerType(String cliName, String description) {
			this.cliName = cliName;
		}

		@Override
		public String toString() { return cliName; }
	}

	static class CallerConverter implements ITypeConverter<CallerType> {
		@Override
		public CallerType convert(String value) throws TypeConversionException {
			String normalized = value.toLowerCase().trim();
			for (CallerType c : CallerType.values()) {
				if (c.toString().equals(normalized)) return c;
			}
			throw new TypeConversionException(
				"Invalid caller '" + value + "'. Valid: auto, ngsep, gatk, freebayes");
		}
	}

	// ── Options ────────────────────────────────────────────────────────────────

	@Option(names = {"-v", "--vcf"},
		description = "Path to the input VCF file.",
		required = true)
	private String vcfFile;

	@Option(names = {"-p", "--ploidy"},
		description = "Ploidy level of the organism (e.g. 2, 4, 10). Default: ${DEFAULT-VALUE}.",
		defaultValue = "2")
	private int ploidy;

	@Option(names = {"-i", "--impute"},
		description = "Imputation method: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}.",
		defaultValue = "bsdp",
		converter = ImputationConverter.class)
	private ImputationMethod impute;

	@Option(names = {"-c", "--caller"},
		description = "Variant caller used to generate the VCF: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}.",
		defaultValue = "auto",
		converter = CallerConverter.class)
	private CallerType caller;

	// ── Execution ──────────────────────────────────────────────────────────────

	@Override
	public Integer call() throws Exception {
		try {
			GeneDosis dosiscgene = new GeneDosis();
			dosiscgene.genDosisAlelicas(vcfFile, ploidy, impute.toInternalKey(), false, caller.toString());
			return 0;
		} catch (Exception e) {
			System.err.println("Error running allele-dosage: " + e.getMessage());
			e.printStackTrace();
			return 1;
		}
	}
}
