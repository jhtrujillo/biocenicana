package huellamolecular;

import java.io.IOException;
import java.text.DecimalFormat;

public class genDosisTargeted {

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
	
	
	
	public void calcularprofunidadTR(String PathFile, int ploidy) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero(PathFile);

		int n = ploidy;
		if (n < 2) {
			n = 2;
		}
		float ploidyLevels[] = new float[n + 1];
		// generate ploidy range for individual from dosage data
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		// ------------------------------------------------------
		// obtengo el nombre de los individuos
		// -------------------------------------------------------
		int cont = 0;
		for (int j = 0; j < ar.numerolineas; j++) {
			String linea = datos[j];
			if (linea.contains("#")) {
				cont = j - 1;
			}
		}
		
		String salida = "CHROM-POS\tREF\tALT";
		for (int i = 9; i < datos[cont + 1].split("	").length; i++) {
			salida = salida + "\t" + datos[cont + 1].split("	")[i];
		}
		System.out.println(salida); // Imprimo la lina de nombres

		// ------------------------------------------------------
		// Obtengo el SNP y el genotipo por cada invidiuo.
		// ------------------------------------------------------
		cont++;
		for (int i = cont+1; i < ar.numerolineas; i++) {

			try {
				String linea = datos[i];
				String[] genotipos = linea.split("	");
				
				String snp = genotipos[0];
				
				String ref = genotipos[3];
				String Alt = genotipos[4];
				String chr = snp.split(":")[0];
				int pos = Integer.parseInt(snp.split(":")[1].split("-")[0]) + 200;
				salida = chr + "_" + pos + "\t" + ref + "\t" + Alt;
				
				
				for (int j = 9; j < genotipos.length; j++) {
					String countRef=genotipos[j].split(":")[3].replace(".","-1");
					String countAlt=genotipos[j].split(":")[5].replace(".","-1");;
					String sumalelos=genotipos[j].split(":")[4].replace(".","-1");;
					String DP = genotipos[j].split(":")[1];
					
					
					salida=salida+"\t"+DP;
				}
				
				
				System.out.println(salida);

			} catch (Exception e) {
				//System.out.println(e);
			}
			
			
			
			
		}

	}
	
	

	public void generarDosis(String PathFile, int ploidy) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero(PathFile);

		int n = ploidy;
		if (n < 2) {
			n = 2;
		}
		float ploidyLevels[] = new float[n + 1];
		// generate ploidy range for individual from dosage data
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		// ------------------------------------------------------
		// obtengo el nombre de los individuos
		// -------------------------------------------------------
		int cont = 0;
		for (int j = 0; j < ar.numerolineas; j++) {
			String linea = datos[j];
			if (linea.contains("#")) {
				cont = j - 1;
			}
		}
		
		String salida = "CHROM-POS\tREF\tALT";
		for (int i = 9; i < datos[cont + 1].split("	").length; i++) {
			salida = salida + "\t" + datos[cont + 1].split("	")[i];
		}
		System.out.println(salida); // Imprimo la lina de nombres

		// ------------------------------------------------------
		// Obtengo el SNP y el genotipo por cada invidiuo.
		// ------------------------------------------------------
		cont++;
		for (int i = cont+1; i < ar.numerolineas; i++) {

			try {
				String linea = datos[i];
				String[] genotipos = linea.split("	");
				
				String snp = genotipos[0];
				String ref = genotipos[3];
				String Alt = genotipos[4];
				String chr = snp.split(":")[0];
				int pos = Integer.parseInt(snp.split(":")[1].split("-")[0]) + 200;
				salida = chr + "_" + pos + "\t" + ref + "\t" + Alt;
				
				
				for (int j = 9; j < genotipos.length; j++) {
					String countRef=genotipos[j].split(":")[3].replace(".","-1");
					String countAlt=genotipos[j].split(":")[5].replace(".","-1");;
					String sumalelos=genotipos[j].split(":")[4].replace(".","-1");;
					
					float dosis=Float.parseFloat(countRef)/(Float.parseFloat(countRef)+Float.parseFloat(countAlt));				
					
					dosis = roundToArray(dosis, ploidyLevels);
					dosis = Float.valueOf(df.format(dosis));
					
					salida=salida+"\t"+dosis;
				}
				
				
				System.out.println(salida);

			} catch (Exception e) {
				//System.out.println(e);
			}
			
			
			
			
		}

	}

		//esta no sirve ya
	public void generarDosis2(String PathFile, int ploidy) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero(PathFile);
		int bandera = 0;

		// Aquí lo que hago es usar un vector para normalizar la dosis alelica, con base
		// a la ploidia del genoma
		// Este es el mimo metodo que utiliza ngsep.
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

			String linea = datos[i];

			if (linea.contains("CHROM") == true) {

				String[] individuos = linea.split("	");

				System.out.print("Chr\tPos\t" + "Ref\t" + "Alt\t");

				for (int j = 9; j < individuos.length; j++) {
					System.out.print(individuos[j] + "\t");
				}

				bandera = 1;
			}

			if (bandera == 1) {

				String[] genotipos = linea.split("	");

				String snp = genotipos[0];
				String ref = genotipos[3];
				String Alt = genotipos[4];

				if (Alt.contains(",") == false) {

					String chr = snp.split(":")[0];
					int pos = 0;

					int lentPos = snp.split(":").length;

					if (lentPos > 1) {
						pos = Integer.parseInt(snp.split(":")[1].split("-")[0]) + 200;
					} else {
						chr = snp.split(":")[0];
					}

					// System.out.print(chr+"\t"+pos+ "\t" + ref + "\t" + Alt + "\t");

					/*
					 * for (int j = 9; j < genotipos.length; j++) {
					 * 
					 * int tamGenotipo = genotipos[j].split(":").length;
					 * 
					 * 
					 * 
					 * if (tamGenotipo == 8) {
					 * 
					 * float countRef = 0; float countAlt = 0; float dosage = 0;
					 * 
					 * 
					 * if (genotipos[j].split(":")[3].compareTo(".") == 0 &&
					 * genotipos[j].split(":")[5].compareTo(".") == 0) { dosage=-1; } else {
					 * 
					 * countRef=Float.parseFloat(genotipos[j].split(":")[3].replace(".","0"));
					 * countAlt=Float.parseFloat(genotipos[j].split(":")[5].replace(".","0"));
					 * 
					 * if((countRef + countAlt) > 0){ dosage = countRef / (countRef + countAlt); }
					 * 
					 * }
					 * 
					 * dosage = roundToArray(dosage, ploidyLevels); dosage =
					 * Float.valueOf(df.format(dosage));
					 * 
					 * System.out.print(dosage+"\t"); //System.out.print(dosage+"["+genotipos[j]+"]"
					 * + "\t"); }
					 * 
					 * else if (tamGenotipo == 6) {
					 * 
					 * float countRef = 0; float countAlt = 0; float dosage = 0;
					 * 
					 * 
					 * if (genotipos[j].split(":")[3].compareTo(".") == 0 &&
					 * genotipos[j].split(":")[5].compareTo(".") == 0) { dosage=-1; } else {
					 * 
					 * countRef=Float.parseFloat(genotipos[j].split(":")[3].replace(".","0"));
					 * countAlt=0;
					 * 
					 * if((countRef + countAlt) > 0){ dosage = countRef / (countRef + countAlt); }
					 * 
					 * if((countRef + countAlt) > 0){ dosage = countRef / (countRef + countAlt); }
					 * 
					 * }
					 * 
					 * dosage = roundToArray(dosage, ploidyLevels); dosage =
					 * Float.valueOf(df.format(dosage));
					 * 
					 * System.out.print(dosage+"\t"); System.out.print(dosage + "\t"); }
					 * 
					 * 
					 * }
					 * 
					 */

					System.out.println("");
				}

			}

		}

	}

	public static void main(String[] args) throws IOException {
		genDosisTargeted dosiscgene = new genDosisTargeted();
		//dosiscgene.generarDosis("/home/estuvar4/Downloads/CEN_132101_FreeBayes_SNPs_Raw.vcf", 10); // imputar usando el
																									// bdsp, quedan -1
																									// si ambos alelos
																									// estan en ceros
		
		dosiscgene.calcularprofunidadTR("/home/estuvar4/Downloads/CEN_132101_FreeBayes_SNPs_Raw.vcf", 10);

	}

}
