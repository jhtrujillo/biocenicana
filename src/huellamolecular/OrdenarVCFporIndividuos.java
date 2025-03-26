package huellamolecular;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class OrdenarVCFporIndividuos {
	archivos ar;
	archivos ar2;
	String[] datos;
	String[] individuos;

	public void loadVCF(String vcfFile) throws IOException {
		ar = new archivos();
		datos = ar.leerfichero2(vcfFile);
		// ar.imprimirDatosFichero();
	}

	public int posInd(String Individuo) {
		int pos = 0;
		for (int i = 0; i < ar.numerolineas; i++) {
			if (datos[i].contains("CHROM") == true) {
				String[] individuos = datos[i].split("	");
				// System.out.println(datos[i]);
				for (int j = 0; j < individuos.length; j++) {
					if (individuos[j].compareTo(Individuo) == 0) {
						pos = j;
					}
				}
				i = ar.numerolineas;
			}
		}
		// System.out.println(pos);
		return pos;
	}

	public String datosSNPvsIndi(String SNP, String Pos, String Individuo) {
		String resultado = "./.:0,0,0:0:0:0,0,0,0:0,0";
		for (int i = 0; i < ar.numerolineas; i++) {
			if (datos[i].contains("#") != true) {
				String snptmp = datos[i].split("	")[0];
				String snppos = datos[i].split("	")[1];

				if (snptmp.compareTo(SNP) == 0 && snppos.compareTo(Pos) == 0) {
					resultado = datos[i].split("	")[posInd(Individuo)];
					if (resultado.contains("/") == false) {
						resultado = "./.:0,0,0:0:0:0,0,0,0:0,0";
					}

				}
			}
		}
		 //System.out.println(resultado);
		return resultado;
	}

	public String evaluarlistadoIndividuos(String Chr, String Pos, String listado) {
		// ar2 = new archivos();
		// individuos = ar2.leerfichero2(listado);

		String resultado = "";

		for (int i = 0; i < ar2.numerolineas; i++) {
			if (i == 0) {
				resultado = datosSNPvsIndi(Chr, Pos, individuos[i]);
			} else {
				resultado = resultado + "	" + datosSNPvsIndi(Chr, Pos, individuos[i]);
			}
		}

		// System.out.println(resultado);
		return resultado;
	}

	public void ordenarVCFxListadoIndividuos(String VCFfile, String listado) {

		ar = new archivos();
		datos = ar.leerfichero2(VCFfile);

		ar2 = new archivos();
		individuos = ar2.leerfichero2(listado);

		for (int i = 0; i < ar.numerolineas; i++) {
			if (datos[i].contains("#CHROM") == true) {
				
				String row=this.datos[i].split("	")[0];
				
				for (int k=1; k<9; k++) {
					row=row+"	"+this.datos[i].split("	")[k];
				}
				row=row+"	";
				//String row = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
				

				for (int j = 0; j < ar2.numerolineas; j++) {
					if (j==0) {
						row = row + individuos[j];
					}else {
						row = row +"	"+ individuos[j];
					}
					
				}
				
				System.out.println(row);
				
			} else if (datos[i].contains("##") == true) {
				System.out.println(datos[i]);
			}
			else {
				String infoSNP=this.datos[i].split("	")[0];
				String Chr = this.datos[i].split("	")[0];
				String Pos = this.datos[i].split("	")[1];
				
				
				for (int k=1; k<9; k++) {
					infoSNP=infoSNP+"	"+this.datos[i].split("	")[k];
				}
			
				
				String genotiposOrdenados = evaluarlistadoIndividuos(Chr, Pos, listado);
				
				System.out.println(infoSNP + "	"+genotiposOrdenados);

			}  

		}
	}

	public static void main(String[] args) throws IOException {
		OrdenarVCFporIndividuos sortVCF = new OrdenarVCFporIndividuos();
		// sortVCF.loadVCF("/home/estuvar4/Downloads/huellaSP803280.vcf");
		// sortVCF.datosSNPvsIndi("SCSP803280_000030297", "6938", "9");

		//sortVCF.loadVCF("/home/estuvar4/Downloads/huellaspontaneum.vcf");
		//sortVCF.datosSNPvsIndi("Chr1B", "11162632", "117");
		
		//sortVCF.loadVCF("/home/estuvar4/Downloads/huellaSorgo.vcf");
		//sortVCF.datosSNPvsIndi("Chr02", "23531219", "104");

		// sortVCF.evaluarlistadoIndividuos("Chr1B",
		// "11162632","/home/estuvar4/Downloads/listado220.txt");
		//sortVCF.ordenarVCFxListadoIndividuos("/home/estuvar4/Downloads/huellaspontaneum.vcf",
		//"/home/estuvar4/Downloads/listado_spontaneum.txt");

		//sortVCF.ordenarVCFxListadoIndividuos("/home/estuvar4/Downloads/huellaSorgo.vcf",
		//"/home/estuvar4/Downloads/listado220.txt");

	
		// sortVCF.posInd("102");

	}

}
