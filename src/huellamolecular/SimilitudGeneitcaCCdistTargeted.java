package huellamolecular;

import java.io.IOException;

public class SimilitudGeneitcaCCdistTargeted {

	Double minDosis = 1.0;
	Double maxDosis = 0.0;

	public void formattoVCF(String pathVCF) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(pathVCF);

		System.out.print("##fileformat=VCFv4.2\n"
				+ "##INFO=<ID=CNV,Number=1,Type=Integer,Description=\"Number of samples with CNVs around this variant\">\n"
				+ "##INFO=<ID=TA,Number=1,Type=String,Description=\"Variant annotation based on a gene model\">\n"
				+ "##INFO=<ID=TID,Number=1,Type=String,Description=\"Id of the transcript related to the variant annotation\">\n"
				+ "##INFO=<ID=TGN,Number=1,Type=String,Description=\"Name of the gene related to the variant annotation\">\n"
				+ "##INFO=<ID=TCO,Number=1,Type=Float,Description=\"One based codon position of the start of the variant. The decimal is the codon position\">\n"
				+ "##INFO=<ID=TACH,Number=1,Type=String,Description=\"Description of the aminoacid change produced by a non-synonymous mutation. String encoded as reference aminoacid, position and mutated aminoacid\">\n"
				+ "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples genotyped\">\n"
				+ "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n"
				+ "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles in called genotypes\">\n"
				+ "##INFO=<ID=AFS,Number=R,Type=Integer,Description=\"Allele counts over the population for all alleles, including the reference\">\n"
				+ "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of variant\">\n"
				+ "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
				+ "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
				+ "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">\n"
				+ "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
				+ "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n"
				+ "##FORMAT=<ID=ADP,Number=R,Type=Integer,Description=\"Counts for observed alleles, including the reference allele\">\n"
				+ "##FORMAT=<ID=BSDP,Number=4,Type=Integer,Description=\"Number of base calls (depth) for the 4 nucleotides in called SNVs sorted as A,C,G,T\">\n"
				+ "##FORMAT=<ID=ACN,Number=R,Type=Integer,Description=\"Predicted copy number of each allele taking into account the prediction of number of copies of the region surrounding the variant\">\n");

		System.out.println(datos[215]);

		for (int i = 216; i < ar.numerolineas; i++) {

			String[] fila = datos[i].split("	");
			
			String CHROM = fila[0];
			String POS = fila[1];
			String ID = fila[2];
			String REF = fila[3];
			String ALT = fila[4];
			String QUAL = fila[5];
			String FILTER = fila[6];
			String INFO = "NS=0;AN=0;AFS=0,0;MAF=0"; //fila[7];
			String FORMAT = "GT:PL:GQ:DP:BSDP:ACN"; //fila[8];

			System.out.print(CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER
					+ "\t" + INFO + "\t" + FORMAT);

			for (int j = 9; j < fila.length; j++) {
				String datosgenotipo = fila[j];
				System.out.print("\t"+tranformarGentoNGSEP(REF, ALT, datosgenotipo));
			}

			System.out.println();
			
		}
	}
	
	public String tranformarGentoNGSEP(String ref, String alt, String genotipo) {
		String result="";
		String GT=genotipo.split(":")[0];
		String PL="0,0,0";
		String GQ="255";
		int DP=0;
		int conteoRef=0;
		int conteAlt=0;
		
		
		if (genotipo.split(":")[1].compareTo(".")!=0){
			DP=Integer.parseInt(genotipo.split(":")[1]);
		}
		else {
			DP=0;
		}
		
		if (genotipo.split(":")[3].compareTo(".")!=0){
			conteoRef=Integer.parseInt(genotipo.split(":")[3]);
		}else {
			conteoRef=0;
		}
		
		if (genotipo.split(":")[5].compareTo(".")!=0){
			//System.out.println("Ok "+genotipo);
			conteAlt=Integer.parseInt(genotipo.split(":")[5]);
		}else {
			conteAlt=DP-conteoRef;
		}
		
		
		int A=0;
		int C=0;
		int G=0;
		int T=0;
		
		if (ref.compareTo("A")==0) {
			A=conteoRef;
		}else if (ref.compareTo("C")==0) {
			C=conteoRef;
		}
		else if (ref.compareTo("G")==0) {
			G=conteoRef;
		}
		else if (ref.compareTo("T")==0) {
			T=conteoRef;
		}
		
		if (alt.compareTo("A")==0) {
			A=conteAlt;
		}else if (alt.compareTo("C")==0) {
			C=conteAlt;
		}
		else if (alt.compareTo("G")==0) {
			G=conteAlt;
		}
		else if (alt.compareTo("T")==0) {
			T=conteAlt;
		}
		
		String BSDP=A+","+C+","+G+","+T;
		return GT+":"+PL+":"+GQ+":"+DP+":"+BSDP;
	}
	
	public void distanciasgeneticas(String pathDosisTargeted) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(pathDosisTargeted);

		for (int i = 1; i < ar.numerolineas; i++) {
			String[] fila = datos[i].split("	");

			for (int j = 4; j < fila.length; j++) {
				double dosis = Double.parseDouble(fila[j]);
				if (minDosis > dosis && dosis > 0) {
					minDosis = dosis;
				}
				if (maxDosis < dosis && dosis > 0) {
					maxDosis = dosis;
				}
				// System.out.print(dosis + " ");
			}

			// System.out.println();

		}
		System.out.println("min:" + minDosis);
		System.out.println("max:" + maxDosis);

	}
	
	

}
