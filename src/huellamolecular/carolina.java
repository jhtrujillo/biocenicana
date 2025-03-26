package huellamolecular;

import java.io.IOException;

public class carolina {
	
	public void leerdatos(String file1, String file2, String genoma) {
		archivos ar1 = new archivos();
		String[] datos1= ar1.leerfichero2(file1);
		
		
		archivos ar2 = new archivos();
		String[] datos2= ar2.leerfichero2(file2);
		
		//ar2.imprimirDatosFichero();
		
		
		for (int i = 0; i< ar1.numerolineas; i++) {
			String genomaRef = datos1[i].split("	")[0];
			String chr = datos1[i].split("	")[1].replace("Chr0", "").replace("Chr", "");
			String snp = datos1[i].split("	")[2];
			
			if (genomaRef.compareTo(genoma)==0) {
				//System.out.println(chr+" "+snp);
				
				for (int j = 0; j< ar2.numerolineas; j++) {
					String genomaRef2 = datos2[j].split("	")[1];
					String chr2 = datos2[j].split("	")[5].replace("Chr0", "").replace("Chr", "");
					String snp2 = datos2[j].split("	")[6];
					
					//System.out.println(chr2+" "+snp2);
					
					if (chr.compareTo(chr2)==0 && snp.compareTo(snp2)==0) {
						System.out.println(datos1[i]+" "+datos2[j]);
					}
				}
				
			}	
		}
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException {
		carolina ca = new carolina();
		//ca.leerdatos("/home/estuvar4/Downloads/marcadores_filtrados_Gwas - Sheet2.tsv", 
		//		"/home/estuvar4/Downloads/01_comparison_different_phenotypes_part1_refSorghum_sdb_produc - TABLA_GENERAL_v3_corrected (1).tsv",
		//	"Sorghum");
		
		
		ca.leerdatos("/home/estuvar4/Downloads/marcadores_filtrados_Gwas - Sheet2.tsv", 
				"/home/estuvar4/Downloads/02_comparison_different_phenotypes_part2_refSpontaneum_sdb_produc - TABLA_GENERAL_spontaneum.tsv",
				"Spontaneum");
	}
	
	
	
}
