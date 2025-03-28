package VCFs;

import java.io.IOException;

import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class Archivos {
	
	
	public static void main(String[] args) {
        try {
            // Crea una instancia del lector, pasando el nombre del archivo VCF
            VCFFileReader reader = new VCFFileReader("/Users/estuvar4/Downloads/tmp2.vcf");
            
           
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
	
	
}
