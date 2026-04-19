package org.cenicana.bio;

import java.util.ArrayList;

public class ProtocolNotes {
	String[] dosis_kasp = null;
	String[] dosis_gbs_rad = null;
	String[] dosis_gbs_rad_wgs = null;
	String[] dosis_targeted_resequencing = null;

	int numSNPs_dosis_kasp = 0;
	int numSNPs_dosis_gbs_rad = 0;
	int numSNPs_gbs_rad_wgs = 0;
	int numSNPstargeted_resequencing = 0;

	ProtocolNotes(String file_dosis_kasp, String file_dosis_gbs_rad, String file_dosis_gbs_rad_wgs,
			String file_dosis_targeted_resequencing) {

		FileUtils ar1 = new FileUtils();
		dosis_kasp = ar1.leerfichero2(file_dosis_kasp);
		numSNPs_dosis_kasp = ar1.numerolineas();
		// System.out.println(numSNPs_dosis_kasp);

		FileUtils ar2 = new FileUtils();
		dosis_gbs_rad = ar2.leerfichero2(file_dosis_gbs_rad);
		numSNPs_dosis_gbs_rad = ar2.numerolineas();
		// System.out.println(numSNPs_dosis_gbs_rad);

		FileUtils ar3 = new FileUtils();
		dosis_gbs_rad_wgs = ar3.leerfichero2(file_dosis_gbs_rad_wgs);
		numSNPs_gbs_rad_wgs = ar3.numerolineas();
		// System.out.println(numSNPs_gbs_rad_wgs);

		FileUtils ar4 = new FileUtils();
		dosis_targeted_resequencing = ar4.leerfichero2(file_dosis_targeted_resequencing);
		numSNPstargeted_resequencing = ar4.numerolineas();
		// System.out.println(numSNPstargeted_resequencing);

		// System.out.println("All data loaded...");

	}

	public ArrayList<String> get_dosis_kasp(String SNP, String ind) {

		ArrayList<String> resultado = new ArrayList<String>();

		for (int i = 0; i < numSNPs_dosis_kasp; i++) {

			String datos[] = dosis_kasp[i].split(",");
			String snp = datos[0];
			String dosis = datos[5];
			String variedad = datos[6];
			String SubjectID = datos[4];
			String indi = datos[6] + "\t" + datos[7] + "\t" + datos[8];

			if (snp.compareTo(SNP) == 0
					&& (ind.compareTo(datos[6]) == 0 || ind.compareTo(datos[7]) == 0 || ind.compareTo(datos[8]) == 0)) {
				resultado.add(SNP + "\t" + indi + "\t" + dosis);
				// System.out.println(SNP + "\t" + SubjectID + "\t" + indi + "\t" + dosis);
			}
		}
		
		if (resultado.size()==0) {
			resultado.add("NA" + "\t" + "NA" + "\t" + "-1");
		}

		return resultado;
	}

	public ArrayList<String> get_dosis_gbs_rad(String SNP, String ind) {

		ArrayList<String> resultado = new ArrayList<String>();

		String[] listadoInd = dosis_targeted_resequencing[0].split("\t");
		int posInd = 0;

		for (int i = 0; i < listadoInd.length; i++) {
			if (listadoInd[i].compareTo(ind) == 0) {
				posInd = i;
				i = listadoInd.length;
			}
		}

		// System.out.println(posInd+" "+listadoInd[posInd]);

		for (int i = 0; i < numSNPs_dosis_gbs_rad; i++) {

			String datos[] = dosis_gbs_rad[i].split("\t");
			String snp = datos[0];
			String dosis = datos[posInd];

			if (snp.compareTo(SNP) == 0) {
				resultado.add(SNP + "\t" + ind + "\t" + dosis);
				// System.out.println(SNP + "\t" + ind + "\t" + "\t" + dosis);
			}

		}
		
		if (resultado.size()==0) {
			resultado.add("NA" + "\t" + "NA" + "\t" + "-1");
		}
		
		return resultado;
	}

	public void get_dosis_kasp_gbs_rad(String SNP, String ind) {

		ArrayList<String> dosis_kasp = get_dosis_kasp(SNP, ind);
		ArrayList<String> dosis_gbs_radseq = get_dosis_gbs_rad(SNP, ind);

		for (int i = 0; i < dosis_kasp.size(); i++) {

			String salida = "";
			salida = dosis_kasp.get(i);

			for (int j = 0; j < dosis_gbs_radseq.size(); j++) {
				salida = salida + "\t" + dosis_gbs_radseq.get(i).split("\t")[1]+ "\t" + dosis_gbs_radseq.get(i).split("\t")[2];
			}

			if (salida.split("\t")[6].compareTo("")!=0 && salida.split("\t")[6].compareTo(SNP)!=0 ) {
				System.out.println(salida);
			}
			salida = "";
		}

	}

	public void get_all_dosis_kasp_gbs_rad() {
		System.out.println("SNP\tsubjectID\tID-Cenicaña\tID-220\tID-220\tDosis-Kasp\tDosis-GBS-Rad");

		for (int i = 1; i < numSNPs_dosis_kasp; i++) {

			try {
				
				String ind = dosis_kasp[i].split(",")[8];
				String snp = dosis_kasp[i].split(",")[0];
				
				if (Integer.parseInt(ind)>0) {
					//System.out.println(snp + " " + ind);
					get_dosis_kasp_gbs_rad(snp, ind);
				}
				
				
				
			} catch (Exception e) {
			}

		}

	}
	
	
	//--------------------------analisis con TR--------------------------
	public ArrayList<String> get_dosis_tr(String SNP, String ind) {

		ArrayList<String> resultado = new ArrayList<String>();

		String[] listadoInd = dosis_targeted_resequencing[0].split("\t");
		int posInd = 0;
		
		//System.out.println(dosis_targeted_resequencing[0]);

		for (int i = 0; i < listadoInd.length; i++) {
			//System.out.println(listadoInd[i]);
			if (listadoInd[i].compareTo(ind) == 0) {
				posInd = i;
				i = listadoInd.length;
			}
		}
		
		//System.out.println(posInd+" "+listadoInd[posInd]);
		 
		for (int i = 0; i < numSNPstargeted_resequencing; i++) {

			String datos[] = dosis_targeted_resequencing[i].split("\t");
			String snp = datos[0];
			String dosis = datos[posInd];

			if (snp.compareTo(SNP) == 0) {
				resultado.add(SNP + "\t" + ind + "\t" + dosis);
				//System.out.println(SNP + "\t" + ind + "\t" + "\t" + dosis);
			}

		}
		
		if (resultado.size()==0) {
			resultado.add("NA" + "\t" + "NA" + "\t" + "-1");
		}
		
		return resultado;
	}
	
	
	public void get_dosis_kasp_tr(String SNP, String ind) {

		ArrayList<String> dosis_kasp = get_dosis_kasp(SNP, ind);
		ArrayList<String> dosis_gbs_radseq = get_dosis_tr(SNP, "CONTROL_"+ind);

		for (int i = 0; i < dosis_kasp.size(); i++) {

			String salida = "";
			salida = dosis_kasp.get(i);

			for (int j = 0; j < dosis_gbs_radseq.size(); j++) {
				salida = salida + "\t" + dosis_gbs_radseq.get(i).split("\t")[1]+ "\t" + dosis_gbs_radseq.get(i).split("\t")[2];
			}
			
			if (salida.split("\t")[6].compareTo("")!=0 && salida.split("\t")[6].compareTo(SNP)!=0 ) {
				System.out.println(salida);
			}
			
			salida = "";
		}

	}
	
	
	public void get_all_dosis_kasp_tr() {
		System.out.println("SNP\tsubjectID\tID-Cenicaña\tID-220\tID-220\tDosis-Kasp\tDosis-GBS-Rad");

		for (int i = 1; i < numSNPs_dosis_kasp; i++) {

			try {
				
				String ind = dosis_kasp[i].split(",")[8];
				String snp = dosis_kasp[i].split(",")[0];
				
				if (Integer.parseInt(ind)>0) {
					//System.out.println(snp + " " + ind);
					get_dosis_kasp_tr(snp, ind);
				}
				
				
				
			} catch (Exception e) {
			}

		}

	}
	
	
	public static void main(String[] args) {

		String file_dosis_kasp = "/home/estuvar4/Downloads/protocol-notes/dosis-intetek.txt";
		String file_dosis_gbs_rad = "/home/estuvar4/Downloads/protocol-notes/AllSamples_genotypes_minMaf_0.01_minC_50_minI_110_bialelicos_dosis_sorgo.txt";
		String file_dosis_gbs_rad_wgs = "/home/estuvar4/Downloads/protocol-notes/cc_01_1940_gbs_radseq_wgs_standarfiltered_d_0_dosisalelicas.txt";
		String file_dosis_targeted_resequencing = "/home/estuvar4/Downloads/protocol-notes/CEN_132101_FreeBayes_SNPs_Raw_dosis.txt";

		ProtocolNotes pn = new ProtocolNotes(file_dosis_kasp, file_dosis_gbs_rad, file_dosis_gbs_rad_wgs,
				file_dosis_targeted_resequencing);

		// pn.get_dosis_kasp_gbs_rad("super_59_423718","103");
		// pn.get_dosis_kasp("super_59_423718", "103");
		// pn.get_dosis_gbs_rad("Chr01_7280", "101");
		//pn.get_dosis_kasp_gbs_rad("super_59_423718", "66");
		//pn.get_all_dosis_kasp_gbs_rad();
		//pn.get_dosis_kasp_tr("Chr06_46251198", "133");
		pn.get_all_dosis_kasp_tr();

	}
}
