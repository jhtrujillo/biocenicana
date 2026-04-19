package org.cenicana.bio;

public class VCF {

	String [] datos = null;
	int numberlines = 0;
	
	
	
	public void loadVCF(String vcfFile) {
		FileUtils ar = new FileUtils();
		datos = ar.leerfichero2(vcfFile);
		this.numberlines = ar.numerolineas();
	}
	
	
	
}
