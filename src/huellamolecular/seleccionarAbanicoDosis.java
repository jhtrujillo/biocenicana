package huellamolecular;

import java.io.IOException;

//Clase que filtrar los snps con base al archivo de las dosis alelicas. 

public class seleccionarAbanicoDosis {

	public boolean continDosis(double dosismin, double dosismax, String[] dosisSnp) {
		boolean resultado = false;
		for (int i = 1; i < dosisSnp.length; i++) {
			if (Double.parseDouble(dosisSnp[i]) >= dosismin && Double.parseDouble(dosisSnp[i]) <= dosismax) {
				resultado = true;
				i = dosisSnp.length;
			}
		}

		return resultado;
	}

	public void loadDosis(String pathfile) {

		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(pathfile);
		int cont = 0;
		int IndexSnp = datos[0].split("	")[0].split("_").length - 1;
						
		
		String snp = "";
		String pos = "";

		for (int i = 1; i < ar.numerolineas; i++) {

			String[] dosisSnp = datos[i].split("	");

			if (continDosis(0.0, 0.2, dosisSnp) && continDosis(0.21, 0.4, dosisSnp) && continDosis(0.41, 0.6, dosisSnp)
					&& continDosis(0.61, 0.8, dosisSnp) && continDosis(0.81, 1.0, dosisSnp))

			{
				
				if (IndexSnp == 0) {
					snp = datos[i].split("	")[0];
					pos= datos[i].split("	")[1];
				} else if (IndexSnp == 1) {
					snp = datos[i].split("	")[0].split("_")[0] + "_" + datos[i].split("	")[0].split("_")[1];
					pos= datos[i].split("	")[1];;
				}
				
								
				System.out.println(snp + "\t" + pos);
				cont++;
			}

		}
		// System.out.println(cont);
	}
	/*
	public static void main(String[] args) throws IOException {

		seleccionarAbanicoDosis abanicodosis = new seleccionarAbanicoDosis();
		abanicodosis.loadDosis("/home/estuvar4/Documents/test/dosis_spon.txt");

	}
*/
}
