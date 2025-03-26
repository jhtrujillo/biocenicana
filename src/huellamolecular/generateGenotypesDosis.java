package huellamolecular;

import java.util.ArrayList;


//Clase que generar las dosis y genotipos en formato txt a partir de un vcf.
public class generateGenotypesDosis {
	public String genenerateNameGenotipes(String header) {
		String names = "CHR_SNP" + "\t";
		String[] genotipos = header.split("\t");
		for (int i = 1; i < genotipos.length; i++) {
			names = names + genotipos[i] + "\t" + "Dosis" + "\t";
			// System.out.println(genotipos[i]+"\t"+"Dosis");
		}
		System.out.println(names);
		return names;
	}

	// Ref = 0, Alt = 1.
	// 0/0 - the sample is homozygous reference
	// 0/1 - the sample is heterozygous, carrying 1 copy of each of the REFand
	// ALT alleles
	// 1/1 - the sample is homozygous alternate
	public String getGenotipo(String ref, String alt, String genotipo) {
		String _genotipo = "(" + genotipo + ")";
		if (genotipo.compareTo("./.") == 0) {
			_genotipo = "./.";
		} else {
			if (genotipo.compareTo("0/1") == 0) {
				_genotipo = ref + "/" + alt;
			} else if (genotipo.compareTo("1/1") == 0) {
				_genotipo = alt + "/" + alt;
			} else if (genotipo.compareTo("0/0") == 0) {
				_genotipo = ref + "/" + ref;
			}
		}
		return _genotipo;
	}

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
			gt = Double.toString(dosis);							
		
		}

		return gt;
	}

	public ArrayList<String> getalelos(String vcfpath) {
		
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfpath);
		ArrayList<String> resultado = new ArrayList<String>();
		String salida = "";
		int splitPosicion = 5;
		//String headGenotypes = datos[18].replace("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ", "#CHROM_SNP\t");
		//genenerateNameGenotipes(headGenotypes);
		
		
		
		for (int i = 19; i < ar.numerolineas; i++) {
			
			
			
			String[] split = datos[i].split("	");
			int numInd = split.length;
		
			String snp = split[1];
			String chr = split[0];
			String ref = split[3];
			String alt = split[4];
			
			try {
				
				
				
				for (int j = 9; j < numInd; j++) {
			
					salida = chr + "_" + snp + "\t";
					
					String genotipo = split[j].split(":")[0];
					if (split[j].split(":").length == 4) {
						splitPosicion = 3;
					} else {
						splitPosicion = 4;
					}
					
					String genotipo_ind_snps = this.getGenotipo(ref, alt, genotipo);
					
					String frecuencia_alelica = getFrecuenciaAlelica(ref, alt, genotipo,
							split[j].split(":")[splitPosicion]);
					
					salida = salida + genotipo_ind_snps + "\t" + frecuencia_alelica + "\t";
					
					//salida = salida + frecuencia_alelica + "\t";
				}
				resultado.add(salida);
				System.out.println(salida);
				salida = "";
			} catch (Exception e) {
			}
		}
		return resultado;
	}
	
	
public ArrayList<String> getAlelovsDosis(String vcfpath) {
		
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfpath);
		ArrayList<String> resultado = new ArrayList<String>();
		String salida1 = "";
		String salida2 = "";
		int splitPosicion = 5;
		//String headGenotypes = datos[18].replace("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ", "#CHROM_SNP\t");
		//genenerateNameGenotipes(headGenotypes);
				
		
		for (int i = 19; i < ar.numerolineas; i++) {
			
			String[] split = datos[i].split("	");
			int numInd = split.length;
		
			String snp = split[1];
			String chr = split[0];
			String ref = split[3];
			String alt = split[4];
			
			try {
				salida1 = chr + "_" + snp + "\t" + ref+ "/" + alt+"\t" ;
				salida2 = chr + "_" + snp + "\t" + ref+ "/" + alt+"\t";
				for (int j = 9; j < numInd; j++) {
					
					String genotipo = split[j].split(":")[0];
					if (split[j].split(":").length == 4) {
						splitPosicion = 3;
					} 
					else if (split[j].split(":").length == 5) {
						splitPosicion = 4;
					}
					else if (split[j].split(":").length == 6) {
						splitPosicion = 4;
					}
					
					String genotipo_ind_snps = this.getGenotipo(ref, alt, genotipo);
					
					
					String frecuencia_alelica = getFrecuenciaAlelica(ref, alt, genotipo,
							split[j].split(":")[splitPosicion]);
					
					salida1 = salida1 + genotipo_ind_snps  + "\t";
					salida2 = salida2 + frecuencia_alelica + "\t";
									
					//salida = salida + frecuencia_alelica + "\t";
				}
				resultado.add(salida1);
				System.out.println(salida1);
				salida1 = "";
				resultado.add(salida2);
				System.out.println(salida2);
				salida2 = "";
			} catch (Exception e) {
			}
		}
		return resultado;
	}

	
	
}