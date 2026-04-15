package org.cenicana.bio;

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

public class AlleleDosageCalculator {

	public VCFFileReader vcfFileReader;
	public Iterator<VCFRecord> iteratorRecords;
	public List<Sample> samples;
	public int genotypePerSamplesComparison[][];
	public int numSNPs = 0;
	public int numGenotypes = 0;
	DecimalFormat df = new DecimalFormat("#.0"); // Dejo las dosis a un solo decimal. 

	String[][] matrizTranspuestaDosis;
	List<String> matrizDosisString = new ArrayList<String>();

	public void computeAlleleDosage(String vcfFile, int ploidy, String impute) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, false, "auto");
	}

	public void computeAlleleDosage(String vcfFile, int ploidy, String impute, boolean storeInMemory) throws IOException {
		computeAlleleDosage(vcfFile, ploidy, impute, storeInMemory, "auto");
	}

	/**
	 * @param callerType Variant caller: "ngsep" | "gatk" | "freebayes" | "auto"
	 */
	public void computeAlleleDosage(String vcfFile, int ploidy, String impute, boolean storeInMemory, String callerType) throws IOException {

		String[] sampleIds = org.cenicana.bio.io.VcfFastReader.getSampleIds(vcfFile);
		numGenotypes = sampleIds.length;

		StringBuilder headerBuilder = new StringBuilder("Chr\tpos\t");
		for (int i = 0; i < sampleIds.length; i++) {
			headerBuilder.append(sampleIds[i]).append("\t");
		}
		String rowHeader = headerBuilder.toString();

		if (storeInMemory) {
			matrizDosisString.add(rowHeader);
		} else {
			System.out.println(rowHeader);
		}

		int n = ploidy;
		if (n < 2) n = 2;
		float ploidyLevels[] = new float[n + 1];
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		Iterable<String[]> blockIterator = org.cenicana.bio.io.VcfFastReader.iterateDataBlocks(vcfFile);

		for (String[] columns : blockIterator) {
			Hashtable<Float, Integer> modaDosisalelicas = new Hashtable<Float, Integer>();
			float promedioDosisalelicas = 0;
			int numIndigenotipado = 0;

			for (int y = 0; y <= n; y++) {
				float dosis = (1.0f / n) * y;
				modaDosisalelicas.put(Float.valueOf(df.format(dosis)), 0);
			}

			String chr = columns[0];
			String pos = columns[1];
			String format = columns.length > 8 ? columns[8] : "";

			String[] formatTokens = format.split(":");
			int adIdx = -1, roIdx = -1, aoIdx = -1, adpIdx = -1, bsdpIdx = -1;
			for (int i = 0; i < formatTokens.length; i++) {
				switch (formatTokens[i]) {
					case "AD":   adIdx   = i; break;
					case "ADP":  adpIdx  = i; break;
					case "RO":   roIdx   = i; break;
					case "AO":   aoIdx   = i; break;
					case "BSDP": bsdpIdx = i; break;
				}
			}

			StringBuilder rowBuilder = new StringBuilder();
			rowBuilder.append(chr).append("\t").append(pos).append("\t");

			float[] parsedDosages = new float[sampleIds.length];
			boolean[] isMissing = new boolean[sampleIds.length];

			int len = Math.min(columns.length - 9, sampleIds.length);
			for (int i = 0; i < len; i++) {
				String genotypeStr = columns[9 + i];
				String[] gtTokens = genotypeStr.split(":");

				float countRef = 0;
				float countAlt = 0;
				boolean foundCounts = false;

				try {
					switch (callerType) {
						case "ngsep":
							// NGSEP: BSDP → [countAlt, countRef, ...]
							if (bsdpIdx != -1 && gtTokens.length > bsdpIdx) {
								String[] bsdp = gtTokens[bsdpIdx].split(",");
								if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
									countAlt = Float.parseFloat(bsdp[0]);
									countRef = Float.parseFloat(bsdp[1]);
									foundCounts = true;
								}
							}
							break;

						case "gatk":
							// GATK: AD → [countRef, countAlt]
							if (adIdx != -1 && gtTokens.length > adIdx) {
								String[] ads = gtTokens[adIdx].split(",");
								if (ads.length >= 2 && !ads[0].equals(".") && !ads[1].equals(".")) {
									countRef = Float.parseFloat(ads[0]);
									countAlt = Float.parseFloat(ads[1]);
									foundCounts = true;
								}
							}
							break;

						case "freebayes":
							// FreeBayes: RO + AO
							if (roIdx != -1 && aoIdx != -1 && gtTokens.length > Math.max(roIdx, aoIdx)) {
								if (!gtTokens[roIdx].equals(".") && !gtTokens[aoIdx].equals(".")) {
									countRef = Float.parseFloat(gtTokens[roIdx]);
									String aoStr = gtTokens[aoIdx];
									if (aoStr.contains(",")) aoStr = aoStr.split(",")[0];
									countAlt = Float.parseFloat(aoStr);
									foundCounts = true;
								}
							}
							break;

						default: // "auto" — prueba en orden de prioridad
							if (bsdpIdx != -1 && gtTokens.length > bsdpIdx) {
								String[] bsdp = gtTokens[bsdpIdx].split(",");
								if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
									countAlt = Float.parseFloat(bsdp[0]);
									countRef = Float.parseFloat(bsdp[1]);
									foundCounts = true;
								}
							} else if (adIdx != -1 && gtTokens.length > adIdx) {
								String[] ads = gtTokens[adIdx].split(",");
								if (ads.length >= 2 && !ads[0].equals(".") && !ads[1].equals(".")) {
									countRef = Float.parseFloat(ads[0]);
									countAlt = Float.parseFloat(ads[1]);
									foundCounts = true;
								}
							} else if (adpIdx != -1 && gtTokens.length > adpIdx) {
								String[] adp = gtTokens[adpIdx].split(",");
								if (adp.length >= 2 && !adp[0].equals(".") && !adp[1].equals(".")) {
									countRef = Float.parseFloat(adp[0]);
									countAlt = Float.parseFloat(adp[1]);
									foundCounts = true;
								}
							} else if (roIdx != -1 && aoIdx != -1 && gtTokens.length > Math.max(roIdx, aoIdx)) {
								if (!gtTokens[roIdx].equals(".") && !gtTokens[aoIdx].equals(".")) {
									countRef = Float.parseFloat(gtTokens[roIdx]);
									String aoStr = gtTokens[aoIdx];
									if (aoStr.contains(",")) aoStr = aoStr.split(",")[0];
									countAlt = Float.parseFloat(aoStr);
									foundCounts = true;
								}
							}
							break;
					}
				} catch (NumberFormatException e) {
					// Dato malformado → queda missing
				}

				float dosage = 0;
				isMissing[i] = true;

				if (foundCounts && (countRef + countAlt) > 0) {
					dosage = countRef / (countRef + countAlt);
					dosage = roundToArray(dosage, ploidyLevels);
					dosage = Float.valueOf(df.format(dosage));

					int conteo = modaDosisalelicas.get(dosage) + 1;
					modaDosisalelicas.put(dosage, conteo);

					parsedDosages[i] = dosage;
					isMissing[i] = false;
					promedioDosisalelicas += dosage;
					numIndigenotipado++;
				}
			}

			float maxvaluedosege = 0;
			float maxkeydosege = 0;
			if (impute.equals("mode") || impute.equals("bsdp-mode")) {
				Set<Float> setOfKeys = modaDosisalelicas.keySet();
				for (Float key : setOfKeys) {
					if (modaDosisalelicas.get(key) >= maxvaluedosege) {
						maxvaluedosege = modaDosisalelicas.get(key);
						maxkeydosege = key;
					}
				}
			}
			float promedio = (numIndigenotipado > 0) ? Float.valueOf(df.format(promedioDosisalelicas / numIndigenotipado)) : 0;

			for (int i = 0; i < len; i++) {
				if (isMissing[i]) {
					if (impute.equals("bsdp")) {
						rowBuilder.append("-1.0\t");
					} else if (impute.equals("mode") || impute.equals("bsdp-mode")) {
						rowBuilder.append(maxkeydosege).append("\t");
					} else if (impute.equals("mean") || impute.equals("bsdp-mean")) {
						rowBuilder.append(promedio).append("\t");
					} else {
						rowBuilder.append("-1.0\t");
					}
				} else {
					rowBuilder.append(parsedDosages[i]).append("\t");
				}
			}

			String row = rowBuilder.toString();

			if (storeInMemory) {
				matrizDosisString.add(row);
			} else {
				System.out.println(row);
			}

			numSNPs++;
		}
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
		AlleleDosageCalculator dosiscgene = new AlleleDosageCalculator();
		dosiscgene.computeAlleleDosage("/Users/estuvar4/Downloads/variants.vcf", 10, "false"); //impute usando el bdsp, quedan -1 si ambos alelos estan
		// en ceros
		// dosiscgene.computeAlleleDosage("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagemode"); //Imputa usando la moda de la dosis
		// dosiscgene.computeAlleleDosage("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosaeaverage"); // imputa usando la mediana - promedio
		// dosiscgene.computeAlleleDosage("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagebsdpmoda"); //imputa usando el bdsp, pero los que quedan -1, los
		// imputa con la moda.
		// dosiscgene.computeAlleleDosage("/home/estuvar4/Downloads/cc-01-1940.vcf", 10,
		// "dosagebsdaverage"); //imputa usando el bdsep, pero los que quedan con -1 los
		// imputa con el promedio.
		//dosiscgene.printDosisMatrix();
		// dosiscgene.TransposeDosisMatrix();
		// dosiscgene.printTransposeDosisMatrix();

		//dosiscgene.get_depth_sequ(
			//	"/home/estuvar4/Downloads/all_220_gbsradwgs_299_parentales_total_519_individuos_22324snps.vcf"); // impute
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