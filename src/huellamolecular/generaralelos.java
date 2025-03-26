package huellamolecular;

import java.io.IOException;
import java.util.ArrayList;

public class generaralelos {

	public String getFrecuenciaAlelica(String ref, String alt, String genotipo, String dsdp) {
		String[] alelosdsdp = dsdp.split(",");
		double refvalue = 0.0;
		double refalt = 0.0;
		String gt = "./.";
		
	
		if (genotipo.compareTo("./.") == 0) {
			gt = "-1";
		} else {
			if (ref.compareTo("A") == 0) {
				refvalue = Double.parseDouble(alelosdsdp[0]);
			} else if (ref.compareTo("C") == 0) {
				refvalue = Double.parseDouble(alelosdsdp[1]);
			} else if (ref.compareTo("G") == 0) {
				refvalue = Double.parseDouble(alelosdsdp[2]);
			} else if (ref.compareTo("T") == 0) {
				refvalue = Double.parseDouble(alelosdsdp[3]);
			}

			if (alt.compareTo("A") == 0) {
				refalt = Double.parseDouble(alelosdsdp[0]);
			} else if (alt.compareTo("C") == 0) {
				refalt = Double.parseDouble(alelosdsdp[1]);
			} else if (alt.compareTo("G") == 0) {
				refalt = Double.parseDouble(alelosdsdp[2]);
			} else if (alt.compareTo("T") == 0) {
				refalt = Double.parseDouble(alelosdsdp[3]);
			}

			double dosis = refvalue / (refvalue + refalt);
			
			//System.out.println(dosis);
			
			if (refvalue==0 && refalt==0){
				gt = "-1";
			}
			else if (dosis >= 0.0 && dosis <= 0.05) {
				gt = "0:1";
			} else if (dosis > 0.05 && dosis <= 0.15) {
				gt = "1:9";
			} else if (dosis > 0.15 && dosis <= 0.25) {
				gt = "2:8";
			} else if (dosis > 0.25 && dosis <= 0.35) {
				gt = "3:7";
			} else if (dosis > 0.35 && dosis <= 0.45) {
				gt = "4:6";
			} else if (dosis > 0.45 && dosis <= 0.55) {
				gt = "1:1";
			} else if (dosis > 0.55 && dosis <= 0.65) {
				gt = "6:4";
			} else if (dosis > 0.65 && dosis <= 0.75) {
				gt = "7:3";
			} else if (dosis > 0.75 && dosis <= 0.85) {
				gt = "8:2";
			} else if (dosis > 0.85 && dosis <= 0.95) {
				gt = "9:1";
			} else if (dosis > 0.95 && dosis <= 1.05) {
				gt = "1:0";
			}
		}

		return gt;
	}
	
	

	public  ArrayList<String> getalelos(String vcfpath, String opc) {
		
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfpath);
		ArrayList<String> resultado = new ArrayList<String>();
		String salida = "";
		int splitPosicion = 5;
		
		boolean escomentario=true;
		int cont=0;
		while(escomentario) {
			if (datos[cont].toString().contains("#")) {
				cont++;
			}
			else {
				escomentario=false;
			}
		}
		
		for (int i = cont; i < ar.numerolineas ; i++) {
			String[] split = datos[i].split("	");
			int numInd = split.length;
			String snp = split[1];
			String chr = split[0];
			String ref = split[3];
			String alt = split[4];
			
			
			
			try {
				salida = snp + "\t" + chr + "\t" + ref + "\t" + alt + "\t";
				for (int j = 9; j < numInd; j++) {
					
					String genotipo = split[j].split(":")[0];
					
					if (split[j].split(":").length==4){
						splitPosicion=3;
					}
					if (split[j].split(":").length==6){
						splitPosicion=4;
					}
					else{
						splitPosicion=5;
					}
					if (opc.compareTo("geno")==0) {
						salida= salida + genotipo+" ";
					}
					else if (opc.compareTo("dosis")==0) {
						salida = salida + getFrecuenciaAlelica(ref, alt, genotipo, split[j].split(":")[splitPosicion]) + "\t" ;
						
					}
				}
				
				resultado.add(salida);
				System.out.println(salida);
				salida = "";
			}
			catch(Exception e){
				System.out.println(e) ; 
			}	
			
		}
		return resultado;
		
	}
	
	

}