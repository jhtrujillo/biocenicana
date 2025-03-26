package huellamolecular;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import ngsep.clustering.DistanceMatrix;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFDistanceMatrixCalculator;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class geneDosisRapidGenomicVCF {
	
	DecimalFormat df = new DecimalFormat("#.0"); // Dejo las dosis a un solo decimal. 
	
	public float roundToArray(float value, float[] array) {

		float best_abs = 1.0f;
		float rounded = value;

		for (int i = 0; i < array.length; i++) {
			if (best_abs > Math.abs(value - array[i])) {
				best_abs = value - array[i];
				rounded = array[i];
			}
		}

		// System.out.println(value+" "+rounded);

		return rounded;
	}


	public void genDosisAlelicas(String vcfFile, int ploidy) throws IOException {
		archivos ar = new archivos();

		String[] datos = ar.leerfichero2(vcfFile);

		// System.out.println(ar.numerolineas);

		// Aquí lo que hago es usar un vector para normalizar la dosis alelica, con base
		// a la ploidia del genoma
		// Este es el mismo metodo que utiliza ngsep.
		int n = ploidy;
		if (n < 2) {
			n = 2;
		}
		float ploidyLevels[] = new float[n + 1];
		// generate ploidy range for individual from dosage data
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		for (int i = 0; i < ar.numerolineas; i++) {

			// Obtengo los nombres de los individuos
			if (datos[i].contains("#CHROM") == true) {
				System.out.println(datos[i].replace("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	", ""));
			}

			// ME salto los comentarios y solo capturo los SNPs.
			if (datos[i].contains("#") != true) {

				String[] row = datos[i].split("	");
				String chr = row[0];
				String snpPos = row[1];
				String ref = row[3];
				String alt = row[4];

				System.out.print(chr + "_" + snpPos + "\t" + ref + "\t" + alt + "\t");

				for (int j = 9; j < row.length; j++) {
					String[] genotipo = row[j].split(":");

					// GT:Genotype
					// DP:Read Depth
					// AD:Number of observation for each allele
					// RO:Reference allele observation count
					// QR:Sum of quality of the reference observations
					// AO:Alternate allele observation count
					// QA:Sum of quality of the alternate observations
					// GL:Genotype Likelihood, log10-scaled likelihoods of the data
					
					//System.out.print(genotipo.length+" "); 

					String GT = genotipo[0];
					String DP = genotipo[1]; // Profunidad de los reads
					String AD = genotipo[2];
					String RO = genotipo[3]; // Conteos del alelo de referencia
					String QR = genotipo[4];
					//String AO = genotipo[5];
					//String QA = genotipo[6];
					//String GL = genotipo[7];

					// Calculo la dosis alelica utilizando la informacion de los conteos de los
					// alelos.
					float dosage = 0;
					if (DP.compareTo(".")!=0 && RO.compareTo(".")!=0) {
						if ((Float.parseFloat(RO) + Float.parseFloat(DP)) > 0) {
							dosage = Float.parseFloat(RO) / Float.parseFloat(DP);
							dosage = roundToArray(dosage, ploidyLevels);
							dosage = Float.valueOf(df.format(dosage));
						} 
						//System.out.print("["+GT+","+dosage+","+RO+","+DP+"]\t");
						System.out.print(dosage+"\t");
					}else {
						//System.out.print("["+GT+","+"-1"+","+RO+","+DP+"]\t");
						System.out.print("-1"+"\t");
					}
						
					
					
					
				}

				System.out.println("");

			}

		}

	}

	public static void main(String[] args) throws IOException {
		geneDosisRapidGenomicVCF dosiscgene = new geneDosisRapidGenomicVCF();
		dosiscgene.genDosisAlelicas(
			"/home/estuvar4/Downloads/tmp.vcf",
			10); // imputar usando el bdsp, quedan -1 si ambos alelos estan

	}

}