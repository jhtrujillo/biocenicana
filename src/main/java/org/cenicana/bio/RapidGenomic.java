package huellamolecular;

import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;

public class RapidGenomic {

	public Hashtable<String, String> cvf = new Hashtable<String, String>();

	public Hashtable<String, String> posicionesSNPs = new Hashtable<String, String>();

	public int leerarchivo(String fastafile, String patron) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(fastafile);

		String tmp = "-" + patron.concat("-");

		String salida = "";
		int posInicial = 0;

		for (int i = 0; i < ar.numerolineas; i++) {
			if (i > 0) {
				salida = salida + datos[i].replace(patron, tmp);
			} else {
				posInicial = Integer.parseInt(datos[i].split(":")[1].split("-")[0]);
				// System.out.println(posInicial);
				// System.out.println(datos[i]);

			}
		}

		int posSNP = posInicial + salida.split("-")[0].length() + patron.length() + 1;
		System.out.print(datos[0] + ":");
		System.out.println(posInicial + salida.split("-")[0].length() + (patron.length() - 1) / 2);

		System.out.println(salida);

		return posSNP;
	}

	//--------------------------------------------------------------------
	// Busca un patron en la cadena del fastafile
	// Al combinar la salida de esta funcion con: | head -n1 | sed "s/:/ /g" | awk '{print $3}' 
	// Me entrega la posición del SNPs en el genoma completo
	//--------------------------------------------------------------------
	public void leerarchivo2(String fastafile, String patron) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(fastafile);

		String tmp = "-" + patron.concat("-");
		
		//System.out.print(tmp);

		String salida = "";
		int posInicial = 0;

		String cadena = "";

		for (int i = 0; i < ar.numerolineas; i++) {
			if (i > 0) {
				//convierto todo el fasta en una sola linea, para poder buscar el patron.
				cadena = cadena + "" + salida + datos[i].replace(patron, tmp);
				//cadena = cadena + "" + salida + datos[i];
				
			} else {
				posInicial = Integer.parseInt(datos[i].split(":")[1].split("-")[0]);
			}
		}


		
		cadena = cadena.replace(patron, tmp);
		
		//System.out.println(cadena);
		
		
		//Aquí hay un posible error, ya que le estaria faltado sumar 1 a la posicion para que quede en todo el SNP
		//Así como esta, esta dando una posicion anterior. Solo una posicion anterior. 
		int posSNP = posInicial + cadena.split("-")[0].length() + patron.length() + 1;
		System.out.print(datos[0] + ":");
		System.out.println(posInicial + cadena.split("-")[0].length() + (patron.length() - 1) / 2);
		System.out.println(cadena);

		
	}

	public void sacarpatron(String fastafile, int ventana) {

		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(fastafile);

		// ar.imprimirDatosFichero();

		String cadena = "";

		for (int i = 0; i < ar.numerolineas; i++) {
			if (i > 0) {
				cadena = cadena + datos[i];
			}
		}

		int PosCenter = ((cadena.length() - 1) / 2) - 1;

		int ini = PosCenter - ventana;
		int fin = PosCenter + ventana;

		System.out.println(PosCenter);

		System.out.println(cadena.substring(PosCenter, PosCenter + 1));

		// for (int i = ini; i < fin; i++) {
		// System.out.print(cadena.substring(i, i+1));
		// }

	}

	// Esta función toma los SNPs que fueron corregidos, por el error del genoma
	// enmarcarado
	// Los SNPs con su posicion corregida esta en posicionesSNPsgenoCompleto
	// El VCF es en donde estan los SNPs con su posicion en el genoma enmascarado
	// Me retorna los SNPs con su posicion corregida en el genoma completo, con sus
	// ref y alt.
	// Este formato es el que nos solicita Rapid Genomics
	public void generarSNPsformatoRG(String posicionesSNPsgenoCompleto, String VCF) {
		FileUtils posiciones = new FileUtils();
		FileUtils vcf = new FileUtils();

		String[] datosposiciones = posiciones.leerfichero2(posicionesSNPsgenoCompleto);
		String[] datosvcf = vcf.leerfichero2(VCF);

		/*
		 * for (int j = 0; j<vcf.numerolineas; j++ ) { if
		 * (datosvcf[j].contains("#")==true) { System.out.println(datosvcf[j]); } else {
		 * j=vcf.numerolineas; } }
		 */

		System.out.println("Pos_Enmas" + "\t" + "Pos_Enmas_VCF" + "\t" + "Pos_Completo" + " " + "Quality\tREF\tALT");

		for (int i = 0; i < posiciones.numerolineas; i++) {
			String[] fila_datosposiciones = datosposiciones[i].split(" ");
			String snpEnmascarado = fila_datosposiciones[0] + "_" + fila_datosposiciones[1];
			String outputposiciones = fila_datosposiciones[0] + "\t" + fila_datosposiciones[4];
			// System.out.println(snpEnmascarado+"
			// "+fila_datosposiciones[0]+"\t"+fila_datosposiciones[4]);

			for (int j = 0; j < vcf.numerolineas; j++) {

				if (datosvcf[j].contains("#") != true) {
					String[] fila_datosvcf = datosvcf[j].split("\t");
					String snpvcf = fila_datosvcf[0] + "_" + fila_datosvcf[1];
					String outputvcf = fila_datosvcf[2] + "\t" + fila_datosvcf[3] + "\t" + fila_datosvcf[4];

					if (snpEnmascarado.compareTo(snpvcf) == 0) {
						System.out.println(snpEnmascarado + "\t" + snpvcf + "\t" + outputposiciones + " " + outputvcf);
						j = vcf.numerolineas;
					}
				}
			}

		}

	}

	// Esta función toma los SNPs que fueron corregidos, por el error del genoma
	// enmarcarado
	// Los SNPs con su posicion corregida esta en posicionesSNPsgenoCompleto
	// El VCF es en donde estan los SNPs con su posicion en el genoma enmascarado
	// Me retorna los SNPs con su posicion corregida en el genoma completo, con sus
	// ref y alt.
	// Este formato es el que nos solicita Rapid Genomics
	public void generarSNPsformatoRG3(String posicionesSNPsgenoCompleto, String VCF) {
		FileUtils posiciones = new FileUtils();
		FileUtils vcf = new FileUtils();

		String[] datosposiciones = posiciones.leerfichero2(posicionesSNPsgenoCompleto);
		String[] datosvcf = vcf.leerfichero2(VCF);
		
		/*
		for (int i = 0; i < posiciones.numerolineas; i++) {
			String[] fila_datosposiciones = datosposiciones[i].split(" ");
			String snpEnmascarado = fila_datosposiciones[0] + "_" + fila_datosposiciones[1];
			String outputposiciones = fila_datosposiciones[0] + " " + fila_datosposiciones[4];
			this.posicionesSNPs.put(snpEnmascarado , outputposiciones);
			//System.out.println(this.posicionesSNPs.get(snpEnmascarado));
		}*/

		System.out.println("Pos_Enmas" + "\t" + "Pos_Enmas_VCF" + "\t" + "Pos_Completo" + " " + "Quality\tREF\tALT");
		
		
		for (int j = 0; j < vcf.numerolineas; j++) {
			if (datosvcf[j].contains("#") != true) {
				String[] fila_datosvcf = datosvcf[j].split("\t");
				String snpvcf = fila_datosvcf[0] + "_" + fila_datosvcf[1];
				String outputvcf = snpvcf+"\t"+fila_datosvcf[2] + "\t" + fila_datosvcf[3] + "\t" + fila_datosvcf[4];
				this.cvf.put(snpvcf , outputvcf);
				//System.out.println(this.cvf.get(snpvcf));
			}
		}
		
		
		
		for (int i = 0; i < posiciones.numerolineas; i++) {
			String[] fila_datosposiciones = datosposiciones[i].split(" ");
			String snpEnmascarado = fila_datosposiciones[0] + "_" + fila_datosposiciones[1];
			String outputposiciones = fila_datosposiciones[0] + "\t" + fila_datosposiciones[4];
			
			String salida=this.cvf.get(snpEnmascarado);
			
			salida = snpEnmascarado+"\t"+outputposiciones+"\t"+salida;
			
			String tmpsalida[] = salida.split("\t");
			
			if (salida.contains("null")==false) {
				//System.out.println(salida);
				System.out.println(tmpsalida[0]+"\t"+tmpsalida[3]+"\t"+tmpsalida[1]+"\t"+tmpsalida[2]+"\t"+tmpsalida[4]+"\t"+tmpsalida[5]+"\t"+tmpsalida[6]);
			}
	
		}
				
	}
	
	
	

	public static void main(String[] args) throws IOException {
		RapidGenomic rg = new RapidGenomic();
		rg.leerarchivo2("/Users/estuvar4/Downloads/Chr05_3731756-3731956.fan","TGCAGGTACTGCTACTCCATC");
		// rg.leerarchivo2("/home/estuvar4/Downloads/test.fasta","TCTTCCCGAGA");
		//rg.generarSNPsformatoRG3("/Users/estuvar4/Downloads/formato_rg/snps_22324_SNPs_corrected_ventana50.txt",

		//		"/Users/estuvar4/Downloads/formato_rg/tmp.vcf");

		

	}

}
