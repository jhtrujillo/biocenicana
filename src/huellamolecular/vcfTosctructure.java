package huellamolecular;

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

public class vcfTosctructure {
	public VCFFileReader vcfFileReader;
	public Iterator<VCFRecord> iteratorRecords;
	public List<Sample> samples;
	public int numGenotypes = 0;
	public int numSNPs = 0;
	String[][] matrixSrtructure = null;

	public void vcfconverTostructure(String vcfFile, int ploidy, String option) throws IOException {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfFile);
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
				String snp = row[0] + "_" + row[1];

				this.matrixSrtructure[0][snpCountpos] = snp;

				for (int snpXgen = 9; snpXgen < row.length; snpXgen++) {
					String genotipo = row[snpXgen];

					String GT = genotipo.split(":")[0];
					String refAlelle = row[3].replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4");
					String altAlelle = row[4].replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4");
					;

					if (GT.compareTo("./.") == 0) {
						refAlelle = altAlelle = "-1";
					}

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

						} else {
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

	public void vcfconverTostructureAlleles(String vcfFile, int ploidy, String option) throws IOException {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfFile);
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
					// System.out.println("S"+row[individuo]+" "+individuo);
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
				String snp = row[0] + "_" + row[1];
		
				

				// Asigno el snp en la primera fila, columna = snpCountpos
				// En snpCountpos llevo el conteo de la columna que le corresponde a cada SNP.
				this.matrixSrtructure[0][snpCountpos] = snp;
				
				//System.out.println(matrixSrtructure[0].length);
				
				//LLeno los valores (genotipo) de las celdas entre individuos vs snp
				for (int snpXgen = 9; snpXgen < row.length; snpXgen++) {
					String genotipo = row[snpXgen];

					String GT = genotipo.split(":")[0];
					String refAlelle = row[3];
					String altAlelle = row[4];
					
					int ref = Integer.parseInt(genotipo.split(":")[5].split(",")[0]);
					int alt = Integer.parseInt(genotipo.split(":")[5].split(",")[1]);
					
					String genotipoalelos="";
					
					for (int ii=0; ii< ref; ii++) {
						genotipoalelos=genotipoalelos+refAlelle;
					}
					
					for (int jj=0; jj< alt; jj++) {
						genotipoalelos=genotipoalelos+altAlelle;
					}
					
					this.matrixSrtructure[snpXgen-8][snpCountpos]=genotipoalelos;
					
				}

				//System.out.println(q"");
				snpCountpos++;
			}
			
			
		}

	}

	public void printMatrix() {
		for (int i = 0; i < this.matrixSrtructure.length; i++) {
			for (int j = 0; j < this.matrixSrtructure[0].length; j++) {
				System.out.print(matrixSrtructure[i][j] + "\t");
			}
			System.out.println();
		}
	}

	public static void main(String[] args) throws IOException {
		vcfTosctructure vcftosctructure = new vcfTosctructure();
		// vcftosctructure.vcfconverTostructure("/home/estuvar4/Downloads/cc-01-1940.vcf",
		// 979,10);
		vcftosctructure.vcfconverTostructureAlleles("/home/estuvar4/Downloads/cc-01-1940.vcf", 10, "ACN");
		vcftosctructure.printMatrix();

	}

}
