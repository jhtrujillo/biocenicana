package huellamolecular;

import java.io.IOException;

public class fingerprint {
	public static void main(String[] args) throws IOException {

		try {
			String opcion = args[0];
			// genera un archivo dosis.txt de la siguiente manera: SCSP803280_000184304_1460
			// 0.578125 0.5466667 0.525 0.50769234
			if (opcion.compareTo("generarDosis") == 0 || opcion.compareTo("1") == 0) {
				try {
					geneDosis dosiscgene = new geneDosis();
					String vcfFile = args[1];
					int ploidy = Integer.parseInt(args[2]);
					String metodoaimputar=args[3];
					dosiscgene.genDosisAlelicas(vcfFile, ploidy, metodoaimputar);
					dosiscgene.printDosisMatrix();
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [generarDosis | 1] [path_vcf] ploidy opciones_imputar:[dosagebsdp|dosagemode|dosaeaverage|dosagebsdpmoda|dosagebsdaverage]> snps_dosis.txt "+e);
				}
			}

			// recibe el archivo dosis.txt generado con generarDosis, y selecciona los SNPs
			// que tienen un abanico en sus dosis alelicas.
			else if (opcion.compareTo("seleccionarDosisAbanico") == 0 | opcion.compareTo("2") == 0) {
				try {
					seleccionarAbanicoDosis abanicodosis = new seleccionarAbanicoDosis();
					abanicodosis.loadDosis(args[1]);
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [seleccionarDosisAbanico|2] [snps_dosis.txt (generado con opción generarDosis)] > snps_dosis_abanico.txt");
				}
			}

			// genera a partir del vcd un txt con los alelos por cada snps: 17120
			// SCSP803280_000007146 C T 1:1 6:4 4:6
			else if (opcion.compareTo("generarAlelosVCF") == 0 | opcion.compareTo("3") == 0) {
				try {
					generaralelos alelos = new generaralelos();
					alelos.getalelos(args[1], "dosis");
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [generarAlelosVCF|3] [path_vcf] > snps_alelos.txt");
				}
			}

			// Recibe el archivos snps_alelos.txt generado por generarAlelosVCF y selecciona
			// los SNPs que cumplen con los filtros de los Australianos.
			if (opcion.compareTo("seleccionarDosisAUS") == 0 || opcion.compareTo("4") == 0) {
				try {
					filtrosalelicosAUS fa = new filtrosalelicosAUS();
					fa.filtrarporClases(args[1], Integer.parseInt(args[2]));
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [seleccionarDosisAUS|4] [snps_alelos.txt (generado con opción generarAlelosVCF)] numIndividuos > snps_dosis_aus.txt");
				}
			}

			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("FiltrarVCF") == 0 || opcion.compareTo("5") == 0) {
				try {
					FiltrarVCF al = new FiltrarVCF();
					al.filtrarvcf(args[1], args[2]);
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [FiltrarVCF|5] [snps_dosis_aus.txt | snps_dosis_abanico.txt] [path_vcf]");
				}
			}
			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("ReducirHuellaVCF") == 0 || opcion.compareTo("6") == 0) {
				try {
					VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
					;
					vcfmatrix.VCFfingerprint(args[1], args[2], 0.0, 0.0, Integer.parseInt(args[3]));
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [ReducirHuellaVCF|6] [path_vcf_original] [path_vcf_filtrado (nombre del archivo resultado)] numSNP (numSnp en huella)");
				}
			}
			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("similitudGeneitcaCCdist") == 0 || opcion.compareTo("7") == 0) {
				try {
					VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
					vcfmatrix.VCFload(args[1]);
					System.out.print("p:10,GD:3 ");
					vcfmatrix.getSimilitudeStats(System.out);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [similitudGeneitcaCCdist | 7] [path_vcf] ");
				}
			}
			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("frecuenciaAlelos") == 0 || opcion.compareTo("8") == 0) {
				try {
					VCFgetHaplotipes vcfcounter = new VCFgetHaplotipes();
					vcfcounter.CounterHaplotipes(args[1]);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [frecuenciaAlelos|8] [path_vcf] ");
				}
			}

			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("genDosisTargeted") == 0 || opcion.compareTo("9") == 0) {
				try {
					genDosisTargeted genDosis = new genDosisTargeted();
					genDosis.generarDosis(args[1]);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [genDosisTargeted|9] [path_vcf] ");
				}
			}

			else if (opcion.compareTo("ComprarDosisHuellavsTargeted") == 0 || opcion.compareTo("10") == 0) {
				try {
					ComprarDosisHuellavsTargeted cdht = new ComprarDosisHuellavsTargeted();
					cdht.comprarindividuos(args[1], args[2]);
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [ComprarDosisHuellavsTargeted|10] dosisSecuenciacion, dosisTargeted, SNPChr, SnpPos");
				}
			} else if (opcion.compareTo("vcf-to-tab-targeted") == 0 || opcion.compareTo("11") == 0) {
				try {
					vcftotabTargeted vcttotab = new vcftotabTargeted();
					vcttotab.vcfToTab(args[1]);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [vcf-to-tab-targeted|11] vcf_targeted");
				}
			} else if (opcion.compareTo("vcftargetedTovcfNGSEP") == 0 || opcion.compareTo("12") == 0) {
				try {
					String dosistargeted = args[1];
					SimilitudGeneitcaCCdistTargeted smgt = new SimilitudGeneitcaCCdistTargeted();
					smgt.formattoVCF(dosistargeted);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [vcftargetedTovcfNGSEP|12] vcf_targeted");
				}
			} else if (opcion.compareTo("printDistanceMatrix") == 0 || opcion.compareTo("13") == 0) {
				try {
					VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
					vcfmatrix.VCFload(args[1]);
					System.out.print("p:10,GD:3 ");
					vcfmatrix.ImprimirMatrix();
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [printDistanceMatrix | 13] [path_vcf] ");
				}
			}
			// genera a partir del vcd un txt con los alelos por cada snps: 17120
			// SCSP803280_000007146 C T 1:1 6:4 4:6
			else if (opcion.compareTo("VCF-to-tab") == 0 | opcion.compareTo("14") == 0) {
				try {
					generaralelos alelos = new generaralelos();
					alelos.getalelos(args[1], "geno");
				} catch (Exception e) {
					System.out.println(
							"Try: java -jar fingerprint.jar [VCF-to-tab|14] [path_vcf] > vcf-to-tab.csv");
				}
			}
			else if (opcion.compareTo("generarDosisTranspuesta") == 0 || opcion.compareTo("15") == 0) {
				try {
					geneDosis dosiscgene = new geneDosis();
					String vcfFile = args[1];
					int ploidy = Integer.parseInt(args[2]);
					String metodoaimputar= args[3];
					dosiscgene.genDosisAlelicas(vcfFile, ploidy, metodoaimputar);
					dosiscgene.TransposeDosisMatrix();
					dosiscgene.printTransposeDosisMatrix();
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [generarDosisTranspuesta | 15 ] [path_vcf] plody > snps_dosis.txt");
				}
			}
			else if (opcion.compareTo("ordenarVCFxListadoIndividuos") == 0 || opcion.compareTo("16") == 0) {
				try {
					OrdenarVCFporIndividuos sortVCF = new OrdenarVCFporIndividuos();
					sortVCF.ordenarVCFxListadoIndividuos(args[1], args[2]);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [ordenarVCFxListadoIndividuos | 16 ] VCFfile ListadoIndividuos");
				}
			}
			else if (opcion.compareTo("vcfToStructure") == 0 || opcion.compareTo("17") == 0) {
				try {
					vcfTosctructure vcftosctructure = new vcfTosctructure();
					int ploidy = Integer.parseInt(args[2]);
					String option = args[3];
					vcftosctructure.vcfconverTostructure(args[1],ploidy, option);
					vcftosctructure.printMatrix();
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [vcfToStructure | 17 ] VCFfile ploidy [ACN]");
				}
			}
			else if (opcion.compareTo("vcfTostructureAlleles") == 0 || opcion.compareTo("18") == 0) {
				try {
					vcfTosctructure vcftosctructure = new vcfTosctructure();
					int ploidy = Integer.parseInt(args[2]);
					String option = args[3];
					String impute = args[4];
					vcftosctructure.vcfconverTostructureAlleles(args[1],ploidy, option, impute);
					vcftosctructure.printMatrixTranspuesta();
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [vcfTostructureAlleles | 18 ] VCFfile ploidy [ACN|dosage] [true|false]");
				}
			}
			else if (opcion.compareTo("addfuntionstogff") == 0 || opcion.compareTo("19") == 0) {
				try {
					addfunctionsgff gff = new addfunctionsgff();
					gff.loadgff(args[1]);
				} catch (Exception e) {
					System.out.println("Try: java -jar fingerprint.jar [addfuntionstogff | 19 ] gffFile");
				}
			}

			// Recibe un listado de snps a seleecionar en el vcf.
			else if (opcion.compareTo("-h") == 0 || opcion.compareTo("-help") == 0 || opcion.compareTo("--help") == 0
					|| opcion.compareTo("") == 0) {
				try {
					System.out.println("Try: java -jar fingerprint.jar [generarDosis | 1] [path_vcf] ploidy opciones_imputar:[dosagebsdp|dosagemode|dosaeaverage|dosagebsdpmoda|dosagebsdaverage]> snps_dosis.txt ");
					System.out.println(
							"Try: java -jar fingerprint.jar [seleccionarDosisAbanico|2] [snps_dosis.txt (generado con opción generarDosis)] > snps_dosis_abanico.txt");
					System.out.println(
							"Try: java -jar fingerprint.jar [generarAlelosVCF|3] [path_vcf] > snps_alelos.txt");
					System.out.println(
							"Try: java -jar fingerprint.jar [seleccionarDosisAUS|4] [snps_alelos.txt (generado con opción generarAlelosVCF)] numIndividuos > snps_dosis_aus.txt");
					System.out.println(
							"Try: java -jar fingerprint.jar [FiltrarVCF|5] [snps_dosis_aus.txt | snps_dosis_abanico.txt] [path_vcf]");
					System.out.println(
							"Try: java -jar fingerprint.jar [ReducirHuellaVCF|6] [path_vcf_original] [path_vcf_filtrado (nombre del archivo resultado)]  numSNP (numSnp en huella)");

					System.out.println("Try: java -jar fingerprint.jar [similitudGeneitcaCCdist|7] [path_vcf] ");

					System.out.println("Try: java -jar fingerprint.jar [frecuenciaAlelos|8] [path_vcf] ");

					System.out.println("Try: java -jar fingerprint.jar [genDosisTargeted|9] [path_vcf] ");

					System.out.println(
							"Try: java -jar fingerprint.jar [ComprarDosisHuellavsTargeted|10] dosisSecuenciacion, dosisTargeted, SNPChr, SnpPos");

					System.out.println("Try: java -jar fingerprint.jar [vcf-to-tab-targeted|11] vcf_targeted");

					System.out.println("Try: java -jar fingerprint.jar [vcftargetedTovcfNGSEP|12] vcf_targeted");

					System.out.println("Try: java -jar fingerprint.jar [printDistanceMatrix | 13] [path_vcf] ");
					
					System.out.println("Try: java -jar fingerprint.jar [VCF-to-tab|14] [path_vcf] > vcf-to-tab.csv");
					
					System.out.println("Try: java -jar fingerprint.jar [ generarDosisTranspuesta | 15 ] [path_vcf] plody > snps_dosis_transpuesta.txt");
					
					System.out.println("Try: java -jar fingerprint.jar [ordenarVCFxListadoIndividuos | 16 ] VCFfile ListadoIndividuos");
					
					System.out.println("Try: java -jar fingerprint.jar [vcfToStructure | 17 ] VCFfile  ploidy [ACN/DSDP/ALL] ");
					
					System.out.println("Try: java -jar fingerprint.jar [vcfTostructureAlleles | 18 ] VCFfile ploidy [ACN|dosage] [true|false]");
					
					System.out.println("Try: java -jar fingerprint.jar [vcfTostructureAlleles | 18 ] VCFfile ploidy [ACN|dosage] [true|false]");
					
					System.out.println("Try: java -jar fingerprint.jar [addfuntionstogff | 19 ] gffFile");
					
					
				} catch (Exception e) {
					System.out.println("Try java -jar fingerprint.jar -help");
				}
			}
		} catch (Exception e) {
			System.out.println("Try java -jar fingerprint.jar -help");
		}

	}

}
