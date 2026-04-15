package huellamolecular;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class VCFgetHaplotipes {

	public void CounterHaplotipes(String pathVCF) throws IOException {
		//VCFFileReader vcfFileReader = new VCFFileReader(
			//	"/home/estuvar4/workspace/huellamolecular/files/huellaRanqueada.vcf");
		VCFFileReader vcfFileReader = new VCFFileReader(pathVCF);
		Iterator<VCFRecord> iteratorRecords = vcfFileReader.iterator();
		List<Sample> samples = vcfFileReader.getHeader().getSamples();
		
		System.out.println("SequenceName" + "\t" + "Position"+"\t"+"undecided"+"\t"+"homozygousReference"+"\t"+"heterozygous"+"\t"+"homozygousVariant");
				
		// Iterate over every variant in VCF file
		while (iteratorRecords.hasNext()) {
			String AlternativeBase = "";
			String Reference = "";
			String SequenceName = "";
			int Position = 0;
			int undecided = 0;
			int homozygousReference = 0;
			int heterozygous = 0;
			int homozygousVariant = 0;

			VCFRecord vcfRecord = iteratorRecords.next();

			// The variant is a SNV
			if (!(vcfRecord.getVariant() instanceof SNV)) {
				// Process only biallelic SNVs
				continue;
			}

			List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls();
			float genotypes[] = new float[genotypeCalls.size()];
			
			
			for (int i = 0; i < genotypeCalls.size(); i++) {
				CalledGenomicVariant call = genotypeCalls.get(i);
				if (!(call instanceof CalledSNV)) {
					continue;
				}
				CalledSNV snv = (CalledSNV) call;
				SequenceName = snv.getSequenceName();
				Position = snv.getPosition();
				if (snv.getGenotype()==-1){
					undecided++;
				}
				else if (snv.getGenotype()==0){
					homozygousReference++;
				}
				else if (snv.getGenotype()==1){
					heterozygous++;
				}
				else if (snv.getGenotype()==2){
					homozygousVariant++;
				}
				//System.out.print(snv.getGenotype()+"\t");
			}
			System.out.println(SequenceName + "\t" + Position+"\t"+undecided+"\t"+homozygousReference+"\t"+heterozygous+"\t"+homozygousVariant);
		}

		vcfFileReader.close();

	}
	
	/*
	public static void main(String[] args) throws IOException {
		VCFgetHaplotipes vcfcounter = new VCFgetHaplotipes();
		//vcfcounter.CounterHaplotipes(args[0]);
		vcfcounter.CounterHaplotipes("/home/estuvar4/Desktop/huellaRanqueada.vcf");
	}*/
	
}
