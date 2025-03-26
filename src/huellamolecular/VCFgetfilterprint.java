package huellamolecular;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import ngsep.clustering.DistanceMatrix;
import ngsep.vcf.VCFDistanceMatrixCalculator;

public class VCFgetfilterprint {

	private VCFDistanceMatrixCalculator vcfdistancematrix;
	private DistanceMatrix matrix; 
	private double[][] distanceMatrix;
	private double minDistance;
	private double maxDistance;
	private double average;
	private double stdev;

	public void VCFload(String VCFpath) throws IOException {
		vcfdistancematrix = new VCFDistanceMatrixCalculator();
		//vcfdistancematrix.setPloidy(2);
		vcfdistancematrix.setPloidy((byte) 10);
		vcfdistancematrix.setDistanceSource(3);
		matrix = vcfdistancematrix.generateMatrix(VCFpath);
		distanceMatrix = matrix.getDistances();
		
		// matrix.printMatrix(System.out);
	}
	
	public void ImprimirMatrix() {
		matrix.printMatrix(System.out);
	}

	public ArrayList<Number> getSimilitudeStats(PrintStream out) {
		ArrayList result = new ArrayList();
		minDistance = 1;
		maxDistance = 0;
		average = 0;
		stdev = 0;

		for (int i = 0; i < distanceMatrix.length; i++) {
			for (int j = 0; j < distanceMatrix[i].length; j++) {
				if (i != j) {
					if (distanceMatrix[i][j] > maxDistance) {
						maxDistance = distanceMatrix[i][j];
					} 
					if (distanceMatrix[i][j] < minDistance) {
						minDistance = distanceMatrix[i][j];
					}
				}
				average = average + distanceMatrix[i][j];
			}
		}

		average = average / (distanceMatrix.length * distanceMatrix[0].length);

		result.add(1 - minDistance);
		result.add(1 - maxDistance);
		result.add(1- average);

		// Calculo la desviacion standar
		for (int i = 0; i < distanceMatrix.length; i++) {
			for (int j = 0; j < distanceMatrix[i].length; j++) {
				stdev = stdev + Math.pow(average - distanceMatrix[i][j], 2);
			}
		}
		stdev = stdev / (distanceMatrix.length * distanceMatrix[0].length);
		stdev = Math.sqrt(stdev);

		result.add(stdev);

		out.println(result);
		return result;
	}

	// Make a copy from the original VCF without a SNP selected.
	public String VCFduplicatefile(String VCFpath, String VCFpathcopy, int toRemove, double minmaf, double maxmaf) {
		archivos ar = new archivos();
		String SNPSelected = null;

		if (ar.existefichero(VCFpathcopy) == 1) {
			ar.eliminarfichero(VCFpathcopy);
		}
		ar.crearfichero(VCFpathcopy);

		String[] vfcFile = ar.leerfichero2(VCFpath);

		for (int i = 0; i < ar.numerolineas; i++) {

			if (toRemove != i) {

				ar.addlineafichereo(VCFpathcopy, vfcFile[i] + "\n");
			} else {

				//System.out.println(vfcFile[i].split("	")[7].split(";").length);
				int tambloque=vfcFile[i].split("	")[7].split(";").length;
				float maf=-1;
				
				if (tambloque==4) {
					maf = Float.parseFloat(vfcFile[i].split("	")[7].split(";")[3].split("=")[1]);
				}else {
					maf = Float.parseFloat(vfcFile[i].split("	")[7].split(";")[2].split("=")[1]);
				}
				
				
				if (minmaf == 0.0 && maxmaf == 0.0) { // selecciono cualquier
														// snps.
					SNPSelected = vfcFile[i];
					System.out.print(SNPSelected.split("	")[0] + "\t" + SNPSelected.split("	")[1] + "\t");
				} else if (maf < minmaf || maf > maxmaf) {// selecciono solo
															// aquellos snps que
															// esten por fuera
															// de los rangos
															// deseados.
					SNPSelected = vfcFile[i];
					System.out.print(SNPSelected.split("	")[0] + "\t" + SNPSelected.split("	")[1] + "\t");
				}

			}
		}
		
		

		return SNPSelected; 
	}

	public void VCFfingerprint(String VCFpath, String VCFpathcopy, double minmaf, double maxmaf, int minSNP) throws IOException {
		System.out.println("FIltrando vcf: "+VCFpath);
		
		archivos ar = new archivos();
		String VCFtmppath = null;
		boolean loop = true;
		String datos[] = ar.leerfichero2(VCFpath);
		
		
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
		
		
		int numSNps=ar.numerolineas-cont;
		int numIntentos=0;
		//ar.DuplicarFichero(VCFpath, VCFpathcopy);
		
		
		while (loop) {
			
			Random r = new Random();
			int Low = cont;
			int High = ar.numerolineas;

			int toRemove = r.nextInt(High - Low) + Low;

			String snpselected = VCFduplicatefile(VCFpath, VCFpathcopy, toRemove, minmaf, maxmaf);
			 
			
			
			if (snpselected != null) {
				VCFtmppath = VCFpathcopy;
				VCFload(VCFtmppath);
				ArrayList stats = getSimilitudeStats(System.out);
				if (Float.parseFloat(stats.get(0).toString()) == 1.0) {
					// loop = false;
					ar.eliminarfichero(VCFpathcopy);
					numIntentos++;
					System.out.print(" Intentos:"+numIntentos);
				} else {
					numIntentos=0;
					ar.eliminarfichero(VCFpath);
					ar.renameFile(VCFpathcopy, VCFpath);
				}
			}
			
			if (numIntentos>50 || ( ar.numerolineas <= (minSNP+cont)) ) {
				loop=false;
			}
			 
		}
		/*
		 * while (loop) { VCFtmppath = VCFpath; VCFload(VCFtmppath); ArrayList stats =
		 * getSimilitudeStats(System.out);
		 * 
		 * if (Float.parseFloat(stats.get(0).toString()) == 1.0) { loop = false; } else
		 * { Random r = new Random(); int Low = 20; int High = 100; int toRemove =
		 * r.nextInt(High - Low) + Low; VCFduplicatefile(VCFtmppath, VCFpathcopy,
		 * toRemove); ar.eliminarfichero(VCFtmppath); ar.renameFile(VCFpathcopy,
		 * VCFpath); } }
	
		 */
	} 

	/*
	public static void main(String[] args) throws IOException {
		VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
		
		//vcfmatrix.VCFload(args[0]);
		//vcfmatrix.VCFload("/home/estuvar4/workarea/30_07_17/huella_650.vcf");		
		//vcfmatrix.VCFload("/home/estuvar4/workarea/30_07_17/220Samples_CTBE_genotypes_hard_filters_sin_erianthus_dosis_filtradas_650.vcf");		
		
		//System.out.print(args[0]+"\t");
		//vcfmatrix.getSimilitudeStats(System.out);
		//System.out.println("p:2,GD:3");
		//vcfmatrix.VCFfingerprint(args[0], args[1], Double.parseDouble(args[2]), Double.parseDouble(args[3]));

		
		//vcfmatrix.VCFfingerprint(
		//"/home/estuvar4/Desktop/pruebas/AllSamples_GBS+Radseq_spon_genotypes_from_unique_standar_filter.vcf",
		//"/home/estuvar4/Desktop/pruebas/AllSamples_GBS+Radseq_spon_genotypes_from_unique_standar_filter_1.vcf", 0.0, 0.0);

	}*/
	
}