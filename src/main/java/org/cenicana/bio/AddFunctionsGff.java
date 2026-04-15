package huellamolecular;

import java.io.IOException;

public class AddFunctionsGff {

	/*
	 * public static void main(String[] args) { AddFunctionsGff gff = new
	 * addfunctionsgff();
	 * gff.loadgff("/home/estuvar4/Downloads/anotacion_funcion.gff");
	 * 
	 * }
	 */

	public void loadgff(String filegff) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(filegff);

		for (int i = 0; i < ar.numerolineas; i++) {
			String linea = datos[i];

			if (linea.contains("gene")) {

				String[] funcion = datos[i + 1].split(";");
				String note = "";

				for (int j = 0; j < funcion.length; j++) {

					String value = funcion[j];
					// System.out.println(value+" "+value.contains("Note"));
					if (value.contains("Note")) {
						note = funcion[j];
					}
				}

				String[] genLinea = linea.split(";");

				for (int j = 0; j < genLinea.length; j++) {
					String value = genLinea[j];

					if (value.contains("Note")) {
						System.out.println(note);
					} else {
						System.out.print(genLinea[j] + ";");
					}
				}

			} else {
				System.out.println(linea);
			}

		}

	}

	public void fixgffformat(String filegff) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(filegff);

		 //for (int i = 0; i < 2; i++) {
		for (int i = 0; i < ar.numerolineas; i++) {

			if (i == 0) {
				System.out.println("##gff-version 3");
				i++;
			}

			String[] geneAtributoV = datos[i].split("\t");
			
//			for (int k = 0; k < geneAtributoV.length; k++) {
//				System.out.println(k+" : "+geneAtributoV[k]);
//			}
			

			if (geneAtributoV.length > 1) {

				String chr = geneAtributoV[0];
				String tool = geneAtributoV[1];
				String type = geneAtributoV[2];
				String posIni = geneAtributoV[3];
				String posFin = geneAtributoV[4];
				String adicionales = geneAtributoV[5] + "\t" + geneAtributoV[6] + "\t" + geneAtributoV[7];
				String funciones = "";
				String atributo = "";
				
				//System.out.println(funciones);

				if (type.compareTo("gene") == 0) {

					//System.out.println();
					
					if (geneAtributoV.length == 9) {
						
						funciones = geneAtributoV[8].replace("\t", " ");
										
						atributo = chr + "\t" + tool + "\t" + type + "\t" + posIni + "\t" + posFin + "\t" + adicionales
								+ "\t" + funciones.replace("\t", " ");
						
						

					} else if (geneAtributoV.length > 9) {

						funciones = geneAtributoV[8].replace("\t", " ").replace(";" + chr, "");

						String atributo1 = chr + "\t" + tool + "\t" + type + "\t" + posIni + "\t" + posFin + "\t"
								+ adicionales + "\t" + funciones.replace("\t", " ");

						String atributo2 = chr + "\t" + tool + "\t" + "mRNA" + "\t" + posIni + "\t" + posFin + "\t"
								+ adicionales + "\t" + funciones.replace("\t", " ");

						atributo = atributo1 + "\n" + atributo2;
					}

					// Evaluo si tiene un Dbxref en la linea siguiente a la del gen.
					String nextgeneAtributoV = datos[i + 1];
					
					
					

					if (nextgeneAtributoV.split(";")[0].contains("Dbxref")) {

						String Dbxref = nextgeneAtributoV.split(";")[0].replace("\t", " ");
						;
						String atributo1 = atributo + " " + Dbxref;
						String atributo2 = atributo + " " + Dbxref;
						
						
						//System.out.println("***** Si lo tiene");

						atributo = atributo1 + "\n" + atributo2.replace("gene", "mRNA");

						i = i + 1;

					}

					System.out.println(atributo);
					

				} else {
					if (geneAtributoV.length == 9) {
						funciones = geneAtributoV[8].replace("\t", " ");
						atributo = chr + "\t" + tool + "\t" + type + "\t" + posIni + "\t" + posFin + "\t" + adicionales
								+ "\t" + funciones.replace("\t", " ");

					}
					System.out.println(atributo);
				}

			}

		}

	}

	public void filtrargffporTamano(String filegff, String minSize, String maxSize) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(filegff);
		// for (int i = 1; i < 20; i++) {
		
		for (int i = 0; i < ar.numerolineas; i++) {
			
			if (i == 0) {
				System.out.println("##gff-version 3");
				i++;
			}

			String[] geneAtributoV = datos[i].split("\t");

			if (geneAtributoV.length > 1) {
				
				
				String chr = geneAtributoV[0];
				String tool = geneAtributoV[1];
				String type = geneAtributoV[2];
				int posIni = Integer.parseInt(geneAtributoV[3]);
				int posFin = Integer.parseInt(geneAtributoV[4]);
				int delta = posFin - posIni;
				String adicionales = geneAtributoV[5] + "\t" + geneAtributoV[6] + "\t" + geneAtributoV[7];
				String funciones =  geneAtributoV[8];
				String atributo = "";
			
				atributo = chr + "\t" + tool + "\t" + type + "\t" + posIni + "\t" + posFin + "\t" + adicionales
						+ "\t" + funciones;
				
				if (type.compareTo("gene") == 0 && (delta >= Integer.parseInt(minSize) && delta <= Integer.parseInt(maxSize) ) ) {
					//System.out.println();
					System.out.println(atributo);
					int flat=1;
					
					while(flat==1 && i < ar.numerolineas-1) {
						//System.out.println(i+" "+ar.numerolineas);
						i++;
						String tmp_atributos = datos[i];
						if (!tmp_atributos.contains("gene")) {
							System.out.println(tmp_atributos);
						}else {
							i--;
							flat=0;
						}
					}
					
					//System.out.println("Termino");

					
				}

			}
			
			
		}

	}
	
	
	public void ajustargfffunciones(String filegff) {
		FileUtils ar = new FileUtils();
		String[] datos = ar.leerfichero2(filegff);
		
		
		for (int i = 0; i < ar.numerolineas; i++) {
			
			if (i == 0) {
				System.out.println("##gff-version 3");
				i++;
			}

			String[] geneAtributoV = datos[i].split("\t");

			if (geneAtributoV.length > 1) {
				
				
				String chr = geneAtributoV[0];
				String tool = geneAtributoV[1];
				String type = geneAtributoV[2];
				int posIni = Integer.parseInt(geneAtributoV[3]);
				int posFin = Integer.parseInt(geneAtributoV[4]);
				int delta = posFin - posIni;
				String adicionales = geneAtributoV[5] + "\t" + geneAtributoV[6] + "\t" + geneAtributoV[7];
				String funciones =  geneAtributoV[8];
				String atributo = "";
			
				atributo = chr + "\t" + tool + "\t" + type + "\t" + posIni + "\t" + posFin + "\t" + adicionales
						+ "\t" + funciones.replace("Note=Protein of unknown function;", "Note=");
				
				
				
				if (type.compareTo("gene") == 0 ) {
					
					String sgt_geneAtributoV = datos[i+1];
					
					String sgt_note=sgt_geneAtributoV.split("Note=")[1];
					
					System.out.println(atributo+sgt_note);
					System.out.println(datos[i+1]);
					
					i++;
					
				}else {
					System.out.println(atributo);
				}
				
				
			}
			
		}
	}

	public static void main(String[] args) throws IOException {
		AddFunctionsGff gff = new AddFunctionsGff();
		//gff.fixgffformat("/Users/estuvar4/Downloads/tmp.gff3");

		//gff.filtrargffporTamano("/Users/estuvar4/Downloads/tmp.gff3", "1000", "5000");
		
		gff.ajustargfffunciones("/Users/estuvar4/Downloads/tmp.gff3");
		
	}

}
