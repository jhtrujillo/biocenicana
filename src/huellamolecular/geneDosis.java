package huellamolecular;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import ngsep.clustering.DistanceMatrix;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFDistanceMatrixCalculator;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class geneDosis {

	public VCFFileReader vcfFileReader;
	public Iterator<VCFRecord> iteratorRecords;
	public List<Sample> samples;
	public int genotypePerSamplesComparison[][];
	public int numSNPs = 0;
	public int numGenotypes = 0;
	Hashtable<String, String[]> hashtable = new Hashtable<String, String[]>();
	String[][] matrizTranspuestaDosis;
	List<String> matrizDosisString = new ArrayList<String>();

	public void genDosisAlelicas(String vcfFile, int ploidy) throws IOException {
		vcfFileReader = new VCFFileReader(vcfFile);
		iteratorRecords = vcfFileReader.iterator();

		samples = vcfFileReader.getHeader().getSamples();

		numGenotypes = samples.size();

		Iterator<VCFRecord> iteratorRecords_tmp = iteratorRecords;

		String row = "Chr\tpos\t";

		// System.out.print("Chr\tpos\t");

		for (int i = 0; i < vcfFileReader.getHeader().getSampleIds().size(); i++) {
			row = row + vcfFileReader.getHeader().getSampleIds().get(i) + "\t";
			// System.out.print(vcfFileReader.getHeader().getSampleIds().get(i)+"\t");
		}

		matrizDosisString.add(row);

		int n = ploidy;
		if (n < 2) {
			n = 2;
		}
		float ploidyLevels[] = new float[n + 1];
		// generate ploidy range for individual from dosage data
		for (int y = 0; y <= n; y++) {
			ploidyLevels[y] = (1.0f / n) * y;
		}

		// System.out.println(row);

		while (iteratorRecords_tmp.hasNext()) {
			VCFRecord vcfRecord = iteratorRecords.next();
			GenomicVariant var = vcfRecord.getVariant();
			String[] alleles = var.getAlleles();
			List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls();

			int posicion_snp = vcfRecord.getFirst();
			String snp_name = vcfRecord.getSequenceName();

			// System.out.print(snp_name+"_"+posicion_snp+"-Dossage_"+snp_name+"_"+posicion_snp+"\t");

			row = snp_name + "\t" + posicion_snp + "\t";

			// System.out.print(snp_name + "\t" + posicion_snp + "\t");

			for (int i = 0; i < genotypeCalls.size(); i++) {
				CalledGenomicVariant call = genotypeCalls.get(i);
				byte[] idxCalledAlleles = call.getIndexesCalledAlleles();
				VariantCallReport report = call.getCallReport();

				if (report == null)
					continue;
				float countRef = report.getCount(alleles[0]);
				float countAlt = report.getCount(alleles[1]);
				float dosage = 0;

				if ((countRef + countAlt) > 0) {
					dosage = countRef / (countRef + countAlt);
					dosage = roundToArray(dosage, ploidyLevels);
					// dosage = (float) (Math.round(dosage * 10.0) / 10.0);

					DecimalFormat df = new DecimalFormat("#.0");
					dosage = Float.valueOf(df.format(dosage));

				}

				if (idxCalledAlleles.length == 0) {
					row = row + -1 + "\t";
					// System.out.print(-1 + "\t");
				} else if (idxCalledAlleles.length == 1) {
					row = row + dosage + "\t";
					// System.out.print(dosage + "\t");
				} else {
					row = row + dosage + "\t";
					// System.out.print(dosage + "\t");
				}

			}
			matrizDosisString.add(row);
			// System.out.println(row);
			// System.out.println();
			numSNPs++;
		}

		// System.out.println(numSNPs);
	}

	public void printDosisMatrix() {
		for (int i = 0; i < this.matrizDosisString.size(); i++) {
			String row=this.matrizDosisString.get(i);
			
			for (int j = 0; j < row.split("\t").length; j++) {
				if (j==0) {
					System.out.print(row.split("\t")[0]+"_"+row.split("\t")[1]+"\t");
					j++;
				}else {
					System.out.print(row.split("\t")[j]+"\t");
				}
				
			}
			System.out.println("");
		}
	}

	public float roundToArray(float value, float[] array) {

		float best_abs = 1.0f;
		float rounded = value;

		for (int i = 0; i < array.length; i++) {
			if (best_abs > Math.abs(value - array[i])) {
				best_abs = value - array[i];
				rounded = array[i];
			}
		}
		return rounded;
	}

	public void TransposeDosisMatrix() {

		int numRows = this.numGenotypes + 2;
		int numCol = this.numSNPs + 1;

		this.matrizTranspuestaDosis = new String[numRows][numCol];

		// System.out.println("Num Fil "+numRows);
		// System.out.println("Num Col "+numCol);

		for (int i = 0; i < numCol; i++) {

			String[] tmp = matrizDosisString.get(i).split("\t");

			for (int j = 0; j < tmp.length; j++) {
				matrizTranspuestaDosis[j][i] = tmp[j];
				// System.out.println(j+" "+i+" "+tmp[j]);
			}

		}

	}

	public void printTransposeDosisMatrix() {
		for (int i = 0; i < this.matrizTranspuestaDosis.length; i++) {
			for (int j = 0; j < this.matrizTranspuestaDosis[0].length; j++) {
				System.out.print(this.matrizTranspuestaDosis[i][j] + "\t");
			}
			System.out.println(" ");
		}
	}
/*
	public static void main(String[] args) throws IOException {
		geneDosis dosiscgene = new geneDosis();
		dosiscgene.genDosisAlelicas("/home/estuvar4/Downloads/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf",
				10); 
		 dosiscgene.printDosisMatrix(); dosiscgene.TransposeDosisMatrix();
		//dosiscgene.printTransposeDosisMatrix();
	}
*/
}