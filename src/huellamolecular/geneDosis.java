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

public class geneDosis {

	public VCFFileReader vcfFileReader;
	public Iterator<VCFRecord> iteratorRecords;
	public List<Sample> samples;
	public int genotypePerSamplesComparison[][];
	public int numSNPs = 0;
	public int numGenotypes = 0;
	DecimalFormat df = new DecimalFormat("#.0"); // Dejo las dosis a un solo decimal. 

	String[][] matrizTranspuestaDosis;
	List<String> matrizDosisString = new ArrayList<String>();

	public void genDosisAlelicas(String vcfFile, int ploidy, String imputar ) throws IOException {

		vcfFileReader = new VCFFileReader(vcfFile); // Cargo el vcf en memoria.
		iteratorRecords = vcfFileReader.iterator(); // Apunto al primer SNPs del vcf, este puntero se mueve al siguiente
													// SNP por cada iteración.
		samples = vcfFileReader.getHeader().getSamples(); // Obtengo los individuos genotipados
		numGenotypes = samples.size(); // Numer de genotipos.

		Iterator<VCFRecord> iteratorRecords_tmp = iteratorRecords;

		// Capturo los nombres de los individuos genotipados. Los uno en una sola
		// variable y los imprimo como cabecera.
		String row = "Chr\tpos\t";
		for (int i = 0; i < vcfFileReader.getHeader().getSampleIds().size(); i++) {
			row = row + vcfFileReader.getHeader().getSampleIds().get(i) + "\t";
			// System.out.print(vcfFileReader.getHeader().getSampleIds().get(i)+"\t");
		}

		matrizDosisString.add(row);

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

		// Itero sobre cada SNP del VCF.
		while (iteratorRecords_tmp.hasNext()) {

			Hashtable<Float, Integer> modaDosisalelicas = new Hashtable<Float, Integer>(); // Primer elemento: dosis
																							// alelica, segundo: vector
																							// [num de repeticiones de
																							// la dosis, promedio en que
																							// se presenta la dosis,
																							// dosis calculada con dsdp]
			float promedioDosisalelicas = 0; // variable usada para sumar las dosis calculadas, luego la uso para
												// dividirla por el total de individuos. Sacando el promedio de las
												// dosis.
			int numIndigenotipado = 0;

			n = ploidy;
			
			if (n < 2) {
				n = 2;
			}
			// Inicializo las tablas hash con el ragno de dosis, teniendo en cuenta la
			// ploidia.
			for (int y = 0; y <= n; y++) {
				float dosis = (1.0f / n) * y;
				modaDosisalelicas.put(Float.valueOf(df.format(dosis)), 0);
			}

			VCFRecord vcfRecord = iteratorRecords.next(); // Me hubico en el SNP siguiente
			GenomicVariant var = vcfRecord.getVariant(); // obtengo las variantes de ese SNPs.
															// Intenrmanete tiene la informacion general del SNP, COmo
															// por ejemplo los alelos que se tiene.
															// Si es bialelico, tendra solo dos alelos

			String[] alleles = var.getAlleles(); // Obtendo en un vector los alelos de referencia y alternativo.
			List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls(); // cargo los datos de cada SNP vs genotipo.
			int posicion_snp = vcfRecord.getFirst(); // Saco la posicino del SNPs.
			String snp_name = vcfRecord.getSequenceName(); // Saco el cromosoma donde esta el SNP
			row = snp_name + "\t" + posicion_snp + "\t";

			// Me muevo sobre cada invidiuo, y saco toda la información que corresponde al
			// genotipo de ese SNP.
			for (int i = 0; i < genotypeCalls.size(); i++) {

				CalledGenomicVariant call = genotypeCalls.get(i); // Genotipo del SNP individuo i.
				VariantCallReport report = call.getCallReport(); // Saco la informacion general de ese genotipo del SNP
																	// en el individuo i.
				if (report == null)
					continue;

				float countRef = report.getCount(alleles[0]); // Saco el conteo del alelo de referencia
				float countAlt = report.getCount(alleles[1]); // Saco el contedo del alelo alternativo
				float dosage = 0;

				
				//System.out.println(row+" "+countRef+" "+countAlt+" "+dosage);
				
				// Calculo la dosis alelica utilizando la informacion de los conteos de los
				// alelos.
				if ((countRef + countAlt) > 0) {
					
					dosage = countRef / (countRef + countAlt);
					
					
					System.out.println(row+" "+countRef+" "+countAlt+" "+dosage);
					
					
					dosage = roundToArray(dosage, ploidyLevels);
					dosage = Float.valueOf(df.format(dosage));		
					
					
					
				}
				
				


				// Cuento cuantas veces se repite cada dosis alelica. Esto es por si deseo
				// imputar con la moda.
				int conteo = modaDosisalelicas.get(dosage) + 1;
				modaDosisalelicas.put(dosage, conteo);

				byte[] idxCalledAlleles = call.getIndexesCalledAlleles(); // Me ayuda a evluar si hay valores perdidos
																			// en el genotipo.
				if (idxCalledAlleles.length == 0) {

					if (imputar.compareTo("dosagebsdp") == 0 || imputar.compareTo("dosagebsdpmoda") == 0
							|| imputar.compareTo("dosagebsdaverage") == 0) { // Imputo con el bsdp. Solo asigna si ambos
																				// alelos son mayores a cero.

						if (countRef > 0 && countAlt > 0) {
							dosage = countRef / (countRef + countAlt);
							dosage = roundToArray(dosage, ploidyLevels);
							dosage = Float.valueOf(df.format(dosage));
							row = row + dosage + "\t";
						} else {
							row = row + -1 + "\t";
						}

					} else { // Si no voy a imputar con BDSP, asigno el valor de -1 por valor perdido.
						row = row + -1 + "\t";
					}

				} else if (idxCalledAlleles.length == 1) {
					row = row + dosage + "\t";
					promedioDosisalelicas = promedioDosisalelicas + dosage;
					numIndigenotipado++;
				} else {
					row = row + dosage + "\t";
					promedioDosisalelicas = promedioDosisalelicas + dosage;
					numIndigenotipado++;
				}

			}

			/****************************************************************************/
			/************ En caso de imputar con moda o con promedio ********************/
			/****************************************************************************/

			if (imputar.compareTo("dosagemode") == 0 || imputar.compareTo("dosagebsdpmoda") == 0) {
				float maxvaluedosege = 0;
				float maxkeydosege = 0;

				Set<Float> setOfKeys = modaDosisalelicas.keySet();
				for (Float key : setOfKeys) {
					if (modaDosisalelicas.get(key) >= maxvaluedosege) {
						maxvaluedosege = modaDosisalelicas.get(key);
						maxkeydosege = key;
					}
				}
				row = row.replaceAll("-1", Float.toString(maxkeydosege));

			}

			else if (imputar.compareTo("dosaeaverage") == 0 || imputar.compareTo("dosagebsdaverage") == 0) {
				float promedio = Float.valueOf(df.format(promedioDosisalelicas / numIndigenotipado));
				row = row.replaceAll("-1", Float.toString(promedio));
			}

			matrizDosisString.add(row);

			numSNPs++;
		}

		// System.out.println(numSNPs);

	}

	public void get_depth_sequ(String vcfFile) throws IOException {

		vcfFileReader = new VCFFileReader(vcfFile); // Cargo el vcf en memoria.
		iteratorRecords = vcfFileReader.iterator(); // Apunto al primer SNPs del vcf, este puntero se mueve al siguiente
													// SNP por cada iteración.
		samples = vcfFileReader.getHeader().getSamples(); // Obtengo los individuos genotipados
		numGenotypes = samples.size(); // Numer de genotipos.

		Iterator<VCFRecord> iteratorRecords_tmp = iteratorRecords;

		// Capturo los nombres de los individuos genotipados. Los uno en una sola
		// variable y los imprimo como cabecera.
		String row = "Chr\tpos\t";
		for (int i = 0; i < vcfFileReader.getHeader().getSampleIds().size(); i++) {
			row = row + vcfFileReader.getHeader().getSampleIds().get(i) + "\t";
			// System.out.print(vcfFileReader.getHeader().getSampleIds().get(i)+"\t");
		}

		matrizDosisString.add(row);

		// Itero sobre cada SNP del VCF.
		while (iteratorRecords_tmp.hasNext()) {

			Hashtable<Float, Integer> modaDosisalelicas = new Hashtable<Float, Integer>(); // Primer elemento: dosis
																							// alelica, segundo: vector
																							// [num de repeticiones de
																							// la dosis, promedio en que
																							// se presenta la dosis,
																							// dosis calculada con dsdp]

			float promedioDosisalelicas = 0; // variable usada para sumar las dosis calculadas, luego la uso para
												// dividirla por el total de individuos. Sacando el promedio de las
												// dosis.
			int numIndigenotipado = 0;

			VCFRecord vcfRecord = iteratorRecords.next(); // Me hubico en el SNP siguiente
			GenomicVariant var = vcfRecord.getVariant(); // obtengo las variantes de ese SNPs.
															// Intenrmanete tiene la informacion general del SNP, COmo
															// por ejemplo los alelos que se tiene.
															// Si es bialelico, tendra solo dos alelos

			String[] alleles = var.getAlleles(); // Obtendo en un vector los alelos de referencia y alternativo.
			List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls(); // cargo los datos de cada SNP vs genotipo.
			int posicion_snp = vcfRecord.getFirst(); // Saco la posicino del SNPs.
			String snp_name = vcfRecord.getSequenceName(); // Saco el cromosoma donde esta el SNP
			row = snp_name + "\t" + posicion_snp + "\t";

			// Me muevo sobre cada invidiuo, y saco toda la información que corresponde al
			// genotipo de ese SNP.
			


			
			
			for (int i = 0; i < genotypeCalls.size(); i++) {
				CalledGenomicVariant call = genotypeCalls.get(i); // Genotipo del SNP individuo i.
				VariantCallReport report = call.getCallReport(); // Saco la informacion general de ese genotipo del SNP
			
				
				//System.out.println(row+" "+alleles.length);
				
				
				
				if (report == null)
					continue;

				float countRef = report.getCount(alleles[0]); // Saco el conteo del alelo de referencia
				float countAlt = report.getCount(alleles[1]); // Saco el contedo del alelo alternativo
				
				
				float total = countRef+countAlt;
				
				row = row + total + "\t";

			}

			
			
			matrizDosisString.add(row);

			numSNPs++;
		}

		// System.out.println(numSNPs);

	}

	public void printDosisMatrix() {
		for (int i = 0; i < this.matrizDosisString.size(); i++) {
			String row = this.matrizDosisString.get(i);

			for (int j = 0; j < row.split("\t").length; j++) {
				if (j == 0) {
					System.out.print(row.split("\t")[0] + "_" + row.split("\t")[1] + "\t");
					j++;
				} else {
					System.out.print(row.split("\t")[j] + "\t");
				}

			}
			System.out.println("");
		}
	}

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

	public void TransposeDosisMatrix() {

		int numRows = this.numGenotypes + 2;
		int numCol = this.numSNPs + 1;

		this.matrizTranspuestaDosis = new String[numRows][numCol];

		// System.out.println("Num Fil "+numRows);
		// System.out.println("Num Col "+numCol);

		for (int i = 0; i < numCol; i++) {

			String[] tmp = matrizDosisString.get(i).split("\t");

			for (int j = 0; j < tmp.length; j++) {
				matrizTranspuestaDosis[j][i] = tmp[j];
				// System.out.println(j+" "+i+" "+tmp[j]);
			}

		}

	}

	public void printTransposeDosisMatrix() {
		for (int i = 0; i < this.matrizTranspuestaDosis.length; i++) {
			for (int j = 0; j < this.matrizTranspuestaDosis[0].length; j++) {
				System.out.print(this.matrizTranspuestaDosis[i][j] + "\t");
			}
			System.out.println(" ");
		}
	}

	public static void main(String[] args) throws IOException {
		geneDosis dosiscgene = new geneDosis();
		dosiscgene.genDosisAlelicas("/Users/estuvar4/Downloads/variants.vcf", 10, "false"); //imputar usando el bdsp, quedan -1 si ambos alelos estan
		// en ceros
		// dosiscgene.genDosisAlelicas("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagemode"); //Imputa usando la moda de la dosis
		// dosiscgene.genDosisAlelicas("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosaeaverage"); // imputa usando la mediana - promedio
		// dosiscgene.genDosisAlelicas("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagebsdpmoda"); //imputa usando el bdsp, pero los que quedan -1, los
		// imputa con la moda.
		// dosiscgene.genDosisAlelicas("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagebsdaverage"); //imputa usando el bdsep, pero los que quedan con -1 los
		// imputa con el promedio.
		//dosiscgene.printDosisMatrix();
		// dosiscgene.TransposeDosisMatrix();
		// dosiscgene.printTransposeDosisMatrix();

		//dosiscgene.get_depth_sequ(
			//	"/home/estuvar4/Downloads/all_220_gbsradwgs_299_parentales_total_519_individuos_22324snps.vcf"); // imputar
																															// usando
																															// el
																															// bdsp,
																															// quedan
																															// -1
																															// si
																															// ambos
																															// alelos
																															// estan
																															// en
																															// ceros
		dosiscgene.printDosisMatrix();

	}

}