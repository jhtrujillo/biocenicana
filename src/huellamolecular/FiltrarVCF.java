package huellamolecular;

import java.util.Hashtable;

public class FiltrarVCF {
	
	public int contInicio=0;

	public Hashtable<String, String> vcfHash = new Hashtable<String, String>();

	public void filtrarvcf(String path_dosis, String vcf) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(path_dosis);

		cargarvcf(vcf);
		
		for (int i = 0; i < this.contInicio; i++) {
			System.out.println(this.vcfHash.get(Integer.toString(i)));
		}

		for (int i = 0; i < ar.numerolineas; i++) {
			String cromosoma = datos[i].split("\t")[0];
			String snp = datos[i].split("\t")[1];		
			//System.out.println(cromosoma+" "+snp+" : "+this.vcfHash.get(cromosoma+"_"+snp));
			String output=this.vcfHash.get(cromosoma+"_"+snp);
			
			//if (output.compareTo("null")!=0) {
				System.out.println(this.vcfHash.get(cromosoma+"_"+snp));
			//}
			
		}
		
	}

	public void cargarvcf(String vcf) {

		vcfHash = new Hashtable<String, String>();

		archivos ar2 = new archivos();
		String[] datosvcf = ar2.leerfichero2(vcf);
		
		boolean escomentario=true;
		int cont=0;
		while(escomentario) {
			if (datosvcf[cont].toString().contains("#")) {
				cont++;
			}
			else {
				escomentario=false;
			}
		}
		
		this.contInicio=cont;
		
		for (int i = 0; i < cont; i++) {
			//System.out.println(datosvcf[i]);
			this.vcfHash.put(Integer.toString(i), datosvcf[i]);
		}

		for (int i = cont; i < ar2.numerolineas; i++) {
			String cromosoma = datosvcf[i].split("	")[0];
			String snp = datosvcf[i].split("	")[1];
			//System.out.println(cromosoma+"_"+snp);
			this.vcfHash.put(cromosoma+"_"+snp , datosvcf[i]);
		}
	}
	

	public static void main(String[] args) {
		FiltrarVCF al = new FiltrarVCF();
		//al.filtrarvcf(args[0], args[1]);

		al.filtrarvcf("/home/estuvar4/Downloads/snps_seleccionar.txt" , "/home/estuvar4/Downloads/cc-01-1940.vcf");

	}

	
}