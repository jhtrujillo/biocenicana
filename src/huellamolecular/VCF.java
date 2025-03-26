package huellamolecular;

public class VCF {

	String [] datos = null;
	int numberlines = 0;
	
	
	
	public void loadVCF(String vcfFile) {
		archivos ar = new archivos();
		datos = ar.leerfichero2(vcfFile);
		this.numberlines = ar.numerolineas();
	}
	
	
	
}
