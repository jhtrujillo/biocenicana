package VCFs;

import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class VCFreader {

	public VCFFileReader reader;
	public String REF;
	public String alleles;
	public short QUAL;

	public static String genotypeType(CalledGenomicVariant call) {
		String type = "";

		if (call.isUndecided()) {
			type = "./.";
		} else if (call.isHeterozygous()) {
			type = "0/1";
		} else if (call.isHomozygous()) {
			type = "0/0";
		} else if (call.isHomozygousReference()) {
			type = "1/1";
		}

		return type;

	}

	public void loadVCF(String VCFfilename) {
		try {
			this.reader = new VCFFileReader(VCFfilename);
		} catch (IOException e) {
			System.err.println("Error al leer el archivo VCF: " + e.getMessage());
		}
	}

	public static void getValuesVariant(VCFRecord record) {

		GenomicVariant variant = record.getVariant();

		String REF = variant.getReference();
		String alleles = Arrays.toString(variant.getAlleles());
		short QUAL = variant.getVariantQS();


		System.out.println(REF+"\t"+alleles+"\t"+QUAL+"\t");

		List<CalledGenomicVariant> calls = record.getCalls();
		for (CalledGenomicVariant call : calls) {

			String sampleId = call.getSampleId();

			String GT = genotypeType(call);
			String Alleles = Arrays.toString(call.getCalledAlleles());
			int gq = call.getGenotypeQuality();
			int dp = call.getTotalReadDepth();
			String bsdps = Arrays.toString(call.getAllCounts());
			String acn = Arrays.toString(call.getAllelesCopyNumber());

			System.out.print(sampleId);
			System.out.print("\t" + GT);
			System.out.print("\t" + Alleles);
			System.out.print("\t" + gq);
			System.out.print("\t" + dp);
			System.out.print("\t" + bsdps);
			System.out.print("\t" + acn);

			System.out.println();

		}

	}

	public static void main(String[] args) {

		String filename = "/Users/estuvar4/Downloads/tmp2.vcf";

		VCFreader vcf = new VCFreader();
		vcf.loadVCF(filename);

		try (VCFFileReader reader = new VCFFileReader(filename)) {

			for (VCFRecord record : reader) {

				getValuesVariant(record);

			}
		} catch (IOException e) {
			System.err.println("Error al leer el archivo VCF: " + e.getMessage());
		}
	}

}
