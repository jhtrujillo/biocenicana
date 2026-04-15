package org.cenicana.bio;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class VcfToStructure2 {
	public VCFFileReader vcfFileReader;
	public Iterator<VCFRecord> iteratorRecords;
	public List<Sample> samples;
	public int numGenotypes = 0;
	public int numSNPs = 0;
	String[][] matrixSrtructure = null;

	public void printMatrix() {
		for (int i = 0; i < this.matrixSrtructure.length; i++) {
			for (int j = 0; j < this.matrixSrtructure[0].length; j++) {
				System.out.print(matrixSrtructure[i][j] + "\t");
			}
			System.out.println();
		}
	}
 
	
	/**
	 * Assign real dosage value depending of ploidy.
	 * 
	 * @param value Dosage ratio calculate from raw base pair count.
	 * @param array Ploidy levels array to fit value.
	 * @throws IOException
	 */
	public float roundToArray(float value, float[] array) {

		float best_abs = 1.0f;
		float rounded = value;

		for (int i = 0; i < array.length; i++) {
			if (best_abs > Math.abs(value - array[i])) {
				best_abs = value - array[i];
				rounded = array[i];
			}
		}
		return rounded;
	}

	
	public void printMatrixTranspuesta() {
		for (int i = 0; i < this.matrixSrtructure[0].length; i++) {
			for (int j = 0; j < this.matrixSrtructure.length; j++) {
				
				if (i==0 && j==0) {
					System.out.print("Marker");
					//System.out.print(matrixSrtructure[j][i]);
				}
				else if (i>0 && j==0) {
					System.out.print(matrixSrtructure[j][i]);
				}
				else {
					System.out.print(";"+matrixSrtructure[j][i]);
				}
				
			}
			System.out.println();
		}	
	}
	
	
	
	public void vcfconverTostructure(String vcfFile, int ploidy, String option) throws IOException {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(vcfFile);
		
		ar.imprimirDatosFichero();
		
		//System.out.println(vcfFile+" "+ar.numerolineas);
		
		
		this.numGenotypes = datos[ar.numerolineas - 1].split("\t").length - 9;
		this.numSNPs = ar.numerolineas;
		int snpCountpos = 1;
		
		
		/*-----------------------------------------------------------------------------
		 * Cuento cuantos comentarios hay en el vcf. 
		 * Debo quitar estos valores, para calcular el total de SNPs. 
		 ----------------------------------------------------------------------------- */
		for (int i = 0; i < ar.numerolineas; i++) {
			if (datos[i].contains("#") == true) {
				this.numSNPs--;
			} else {
				i = ar.numerolineas;
			}
		}
		
		
		matrixSrtructure = new String[(this.numGenotypes * ploidy) + 1][this.numSNPs + 1];
		this.matrixSrtructure[0][0] = "";
		
		
		/*-----------------------------------------------------------------------------
		 * Recorro cada fila del VCF. Salto los comentarios
		 * 1. con CHROM saco los individuos
		 * 2. Todo lo que no tenga # es SNPs. 
		 ----------------------------------------------------------------------------- */
		for (int i = 0; i < ar.numerolineas; i++) {

			/*-----------------------------------------------------------------------------
			 * Saco los individuos, los asigno a cada fila en la primera columna. 
			 * Tengo en cuenta la ploidia del genoma
			 ----------------------------------------------------------------------------- */
			if (datos[i].contains("#CHROM") == true) {
				String[] row = datos[i].split("\t"); // Esta linea contiene los individuos
				int posRow = 1;
				for (int individuo = 9; individuo < row.length; individuo++) {
					for (int k = posRow; k < (posRow + ploidy); k++) {
						// System.out.println("S"+row[individuo]+" "+individuo);
						this.matrixSrtructure[k][0] = "S" + row[individuo];
					}
					posRow = posRow + ploidy;
				}
				// System.out.println();
			}
			
			/*-----------------------------------------------------------------------------
			 * Saco los SNPs y los pongo como columnas. Todo lo que no tenga # es SNP.
			 * Tengo en cuenta la ploidia del genoma
			 ----------------------------------------------------------------------------- */
			int posRow = 1;
			if (datos[i].contains("#") == false) {
				String[] row = datos[i].split("\t"); // cada row contiene el SNPs y su genotipo en los individuos.
				String snp = row[0] + "_" + row[1]+"_"+row[3]+"_"+row[4]; //El SNP quedaría con SNP_ref_alt

				this.matrixSrtructure[0][snpCountpos] = snp;

				for (int snpXgen = 9; snpXgen < row.length; snpXgen++) {
					
					String genotipo = row[snpXgen];

					String GT = genotipo.split(":")[0];
					String refAlelle = row[3].replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4");
					String altAlelle = row[4].replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4");
					;

					//if (GT.compareTo("./.") == 0) {
					//	refAlelle = altAlelle = "-1";
					//}

					int cont = 0;
					for (int snpXgenxdosis = posRow; snpXgenxdosis < (posRow + ploidy); snpXgenxdosis++) {

						if (option.compareTo("ACN") == 0) {

							int ref = Integer.parseInt(genotipo.split(":")[5].split(",")[0]);
							int alt = Integer.parseInt(genotipo.split(":")[5].split(",")[1]);

							if (cont < ref) {
								this.matrixSrtructure[snpXgenxdosis][snpCountpos] = refAlelle;
								cont++;
							} else {
								this.matrixSrtructure[snpXgenxdosis][snpCountpos] = altAlelle;
							}

						} 
						else if (option.compareTo("BSDP") == 0) {
							
							String aleloref=row[3];
							String aleloalt=row[4];
							
							int ref = 0;
							int alt = 0;
							
							
							if (aleloref=="A") {
								ref=Integer.parseInt(genotipo.split(":")[4].split(",")[0]);
							}
							if (aleloref=="C") {
								ref=Integer.parseInt(genotipo.split(":")[4].split(",")[1]);
							}
							if (aleloref=="G") {
								ref=Integer.parseInt(genotipo.split(":")[4].split(",")[2]);
							}
							if (aleloref=="T") {
								ref=Integer.parseInt(genotipo.split(":")[4].split(",")[3]);
							}
							
							
							if (aleloalt=="A") {
								alt=Integer.parseInt(genotipo.split(":")[4].split(",")[0]);
							}
							if (aleloalt=="C") {
								alt=Integer.parseInt(genotipo.split(":")[4].split(",")[1]);
							}
							if (aleloalt=="G") {
								alt=Integer.parseInt(genotipo.split(":")[4].split(",")[2]);
							}
							if (aleloalt=="T") {
								alt=Integer.parseInt(genotipo.split(":")[4].split(",")[3]);
							}
							
							
							
							if (cont < ref) {
								this.matrixSrtructure[snpXgenxdosis][snpCountpos] = refAlelle;
								cont++;
							} else {
								this.matrixSrtructure[snpXgenxdosis][snpCountpos] = altAlelle;
							}
							

						}
						else {
							this.matrixSrtructure[snpXgenxdosis][snpCountpos] = genotipo;
						}

					}

					posRow = posRow + ploidy;

				}

				// System.out.println("");
				snpCountpos++;
			}

		}
	}

	
	public void vcfconverTostructureAlleles(String vcfFile, int ploidy, String option, String imputar)
			throws IOException {
		
		
		
		
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(vcfFile);
		
		//System.out.println(ar.numerolineas);
		
		this.numGenotypes = datos[ar.numerolineas - 1].split("\t").length - 9;
		this.numSNPs = ar.numerolineas;
		int snpCountpos = 1;

		
		//System.out.println(numGenotypes);
		
		if (imputar.compareTo("") == 0) {
			imputar = "false";
		}

		int n = ploidy;
		if (n < 2) {
			n = 2;
		}
		float ploidyLevels[] = new float[n + 1];

		// generate ploidy range for individual from dosage data
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
			//System.out.println(ploidyLevels[y]);
		}

		/*-----------------------------------------------------------------------------
		 * Cuento cuantos comentarios hay en el vcf. 
		 * Debo quitar estos valores, para calcular el total de SNPs. 
		 ----------------------------------------------------------------------------- */
		for (int i = 0; i < ar.numerolineas; i++) {
			if (datos[i].contains("#") == true) {
				this.numSNPs--;
			} else {
				i = ar.numerolineas;
			}
		}

		matrixSrtructure = new String[this.numGenotypes + 1][this.numSNPs + 1];
		this.matrixSrtructure[0][0] = "";

		/*-----------------------------------------------------------------------------
		 * Recorro cada fila del VCF. Salto los comentarios
		 * 1. con CHROM saco los individuos
		 * 2. Todo lo que no tenga # es SNPs. 
		 ----------------------------------------------------------------------------- */
		for (int i = 0; i < ar.numerolineas; i++) {

			/*-----------------------------------------------------------------------------
			 * Saco los individuos, los asigno a cada fila en la primera columna. 
			 * Tengo en cuenta la ploidia del genoma
			 ----------------------------------------------------------------------------- */
			if (datos[i].contains("#CHROM") == true) {
				String[] row = datos[i].split("\t"); // Esta linea contiene los individuos
				int posRow = 1;
				for (int individuo = 9; individuo < row.length; individuo++) {
					//System.out.println("S"+row[individuo]+" "+individuo);
					this.matrixSrtructure[individuo - 8][0] = "S" + row[individuo];

				}
				// System.out.println();
			}

			/*-----------------------------------------------------------------------------
			 * Saco los SNPs y los pongo como columnas. Todo lo que no tenga # es SNP.
			 * Proceso un SNPs por cada iteracion con la variable i
			 * Tengo en cuenta la ploidia del genoma
			 ----------------------------------------------------------------------------- */
			int posRow = 1;
			if (datos[i].contains("#") == false) {
				String[] row = datos[i].split("\t"); // cada row contiene el SNPs y su genotipo en los individuos.
				//String snp = row[0] + "_" + row[1];
				String snp = row[0] + "_" + row[1]+"_"+row[3]+"_"+row[4];

				// Asigno el snp en la primera fila, columna = snpCountpos
				// En snpCountpos llevo el conteo de la columna que le corresponde a cada SNP.
				this.matrixSrtructure[0][snpCountpos] = snp;

				//System.out.println(row.length);

				// LLeno los valores (genotipo) de las celdas entre individuos vs snp
				for (int snpXgen = 9; snpXgen < row.length; snpXgen++) {
					String genotipo = row[snpXgen];

					if (option.compareTo("ACN") == 0) {
						String GT = genotipo.split(":")[0];
						String refAlelle = row[3];
						String altAlelle = row[4];

						int ref = Integer.parseInt(genotipo.split(":")[5].split(",")[0]);
						int alt = Integer.parseInt(genotipo.split(":")[5].split(",")[1]);

						String genotipoalelos = "";

						for (int ii = 0; ii < ref; ii++) {
							genotipoalelos = genotipoalelos + refAlelle;
						}

						for (int jj = 0; jj < alt; jj++) {
							genotipoalelos = genotipoalelos + altAlelle;
						}

						if (GT.compareTo("./.") == 0 && imputar.compareTo("false") == 0) {
							genotipoalelos = "NA";
						}

						this.matrixSrtructure[snpXgen - 8][snpCountpos] = genotipoalelos;
						
					} 
					
					else if (option.compareTo("dosage") == 0) {
						String GT = genotipo.split(":")[0];
						String BSDP = genotipo.split(":")[4];// A,C,G,T
						String refAlelle = row[3];
						String altAlelle = row[4];
						float refAlelleCount = 0;
						float altAlelleCount = 0;

						if (refAlelle.compareTo("A") == 0) {
							refAlelleCount = Float.parseFloat(BSDP.split(",")[0]);
						} else if (refAlelle.compareTo("C") == 0) {
							refAlelleCount = Float.parseFloat(BSDP.split(",")[1]);
						} else if (refAlelle.compareTo("G") == 0) {
							refAlelleCount = Float.parseFloat(BSDP.split(",")[2]);
						} else if (refAlelle.compareTo("T") == 0) {
							refAlelleCount = Float.parseFloat(BSDP.split(",")[3]);
						}

						if (altAlelle.compareTo("A") == 0) {
							altAlelleCount = Float.parseFloat(BSDP.split(",")[0]);
						} else if (altAlelle.compareTo("C") == 0) {
							altAlelleCount = Float.parseFloat(BSDP.split(",")[1]);
						} else if (altAlelle.compareTo("G") == 0) {
							altAlelleCount = Float.parseFloat(BSDP.split(",")[2]);
						} else if (altAlelle.compareTo("T") == 0) {
							altAlelleCount = Float.parseFloat(BSDP.split(",")[3]);
						}
						
						float dosage = -1;
						
						if (GT.compareTo("./.")!=0) {
							dosage = refAlelleCount / (refAlelleCount + altAlelleCount);
							dosage = roundToArray(dosage, ploidyLevels);
						}
						
						int ref = 0;
						int alt = 0;

						for (int k = 0; k<ploidyLevels.length; k++) {
							if (dosage==ploidyLevels[k]) {
								ref=k;
								alt=ploidy-ref;
							}
						}
						
						
						String genotipoalelos = "";

						//Consulto si voy a imputar con ACN. false = no voy a imputar
						//if (GT.compareTo("./.") == 0 && imputar.compareTo("false") == 0) {
						//	genotipoalelos = "";
						//}
						//else 
						if (GT.compareTo("./.") == 0 && imputar.compareTo("true") == 0) {
							GT = genotipo.split(":")[0];
							
							//System.out.println(genotipo);
							
							refAlelle = row[3]; 
							altAlelle = row[4];
							
							ref = Integer.parseInt(genotipo.split(":")[5].split(",")[0]); //imputo con el ACN
							alt = Integer.parseInt(genotipo.split(":")[5].split(",")[1]); //imputo con el ACN
							
						}
						 
						
						for (int ii = 0; ii < ref; ii++) {
							genotipoalelos = genotipoalelos + refAlelle;
						}
						
						for (int jj = 0; jj < alt; jj++) {
							genotipoalelos = genotipoalelos + altAlelle;
						}
						
						
						
						this.matrixSrtructure[snpXgen - 8][snpCountpos] = genotipoalelos;
						
						//this.matrixSrtructure[snpXgen - 8][snpCountpos] = GT+" "+Float.toString(dosage)+" "+ref+" "+alt;
					    
					}
				}

				// System.out.println(q"");
				snpCountpos++;
			}

		}

	}
	
	public static void main(String[] args) throws IOException {
		VcfToStructure2 vcftosctructure = new VcfToStructure2();
		//vcftosctructure.vcfconverTostructure("/home/estuvar4/Downloads/tmp.vcf",10, "");
		vcftosctructure.vcfconverTostructureAlleles("/Users/estuvar4/Downloads/tmp.vcf", 10, "dosage", "false");
		vcftosctructure.printMatrixTranspuesta();

	}

}
