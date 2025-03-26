package huellamolecular;

public class vcftotabTargeted {

	public void vcfToTab(String PathFile) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(PathFile);
		int bandera = 0;

		for (int i = 0; i < ar.numerolineas; i++) {
			String linea = datos[i];
			// System.out.println(linea);

			if (linea.contains("CHROM") == true) {

				String[] individuos = linea.split("	");

				System.out.print("Chr\tPos\t" + "Ref\t" + "Alt\t");

				for (int j = 9; j < individuos.length; j++) {
					System.out.print(individuos[j] + "\t");
				}

				bandera = i+1;

			}
		}

		System.out.println("NumNoGenotipados");

		for (int i = bandera; i < ar.numerolineas; i++) {
			String linea = datos[i];

			String[] genotipos = linea.split("	");

			try {
				String snp = genotipos[0];
				String ref = genotipos[3];
				String Alt = genotipos[4];
				String chr = snp.split(":")[0];
				int lentPos = snp.split(":").length;
				int pos = Integer.parseInt(snp.split(":")[1].split("-")[0]) + 200;
				String listadoGenotipos = "";
				int numNoGenotipados = 0;

				System.out.print(chr + "\t" + pos + "\t" + ref + "\t" + Alt + "\t");

				for (int j = 9; j < genotipos.length; j++) {

					int tamGenotipo = genotipos[j].split(":").length;
					String genotipo = genotipos[j].split(":")[0];
					System.out.print(genotipo + "\t");
					if (genotipo.compareTo("./.") == 0) {
						numNoGenotipados++;
					}

				}

				System.out.println(numNoGenotipados + "");

			} catch (Exception e) {
				System.out.println("Error en "+linea);
			}

		}
	}

	/*
	public static void main(String[] args) {
		String vcfFile = "/home/estuvar4/eclipse-workspace/ensamblaje/files/CEN_132101_FreeBayes_SNPs_Raw.vcf";
		vcftotabTargeted vcttotab = new vcftotabTargeted();
		vcttotab.vcfToTab(vcfFile);
	}*/

}
