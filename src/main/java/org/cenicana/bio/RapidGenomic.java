package org.cenicana.bio;

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



	public static void main(String[] args) throws IOException {
		RapidGenomic rg = new RapidGenomic();
		rg.leerarchivo2("/Users/estuvar4/Downloads/Chr05_3731756-3731956.fan","TGCAGGTACTGCTACTCCATC");
	}

}
