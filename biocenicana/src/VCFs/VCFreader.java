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


	public static String genotypeType(CalledGenomicVariant call) {
		String type="";

		if (call.isUndecided()) {
			type="./.";
		}else if (call.isHeterozygous()) {
			type="0/1";
		}
		else if (call.isHomozygous()) {
			type="0/0";
		}
		else if (call.isHomozygousReference()) {
			type="1/1";
		}

		return type;

	}



	public static void main(String[] args) {

		String filename = "/Users/estuvar4/Downloads/tmp2.vcf";
		try (VCFFileReader reader = new VCFFileReader(filename)) {


			for (VCFRecord record : reader) {

				GenomicVariant variant = record.getVariant();
				System.out.println("Variante: " + variant);


				List<CalledGenomicVariant> calls = record.getCalls();
				for (CalledGenomicVariant call : calls) {
					
				    String sampleId = call.getSampleId();
					
					String GT = genotypeType(call);
					String Alleles = Arrays.toString(call.getCalledAlleles());
                    int gq = call.getGenotypeQuality();
                    int dp = call.getTotalReadDepth();
                    String bsdps = Arrays.toString(call.getAllCounts());
                    String acn = Arrays.toString(call.getAllelesCopyNumber());
            
//                    System.out.print("  Muestra: " + sampleId);
//                    System.out.print("    GT: " + GT);
//                    System.out.print("    Alleles: " + Alleles);
//                    System.out.print("    GQ: " + gq);
//                    System.out.print("    DP: " + dp);
//                    System.out.print("    BSDP: " + bsdps);
//                    System.out.print("    ACN: " + acn);
                    
                    
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
		} catch (IOException e) {
			System.err.println("Error al leer el archivo VCF: " + e.getMessage());
		}
	}


}
