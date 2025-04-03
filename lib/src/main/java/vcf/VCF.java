package vcf;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.util.Iterator;
import java.util.Map;



	public class VCF {

	    private VCFHeader header;
	    private Iterator<VariantContext> variantsIterator;

	    // Constructor que lee el archivo VCF
	    public VCF(String vcfFilePath) {
	        File vcfFile = new File(vcfFilePath);
	        VCFFileReader reader = new VCFFileReader(vcfFile, false);

	        this.header = reader.getFileHeader();
	        this.variantsIterator = reader.iterator();
	    }

	    // Método para obtener el encabezado
	    public VCFHeader getHeader() {
	        return header;
	    }

	    // Método para obtener el siguiente VariantContext (una variante)
	    public VariantContext getNextVariant() {
	        if (variantsIterator.hasNext()) {
	            return variantsIterator.next();
	        }
	        return null;  // Si no hay más variantes
	    }

	    // Método para iterar sobre todas las variantes
	    public void iterateVariants() {
	        while (variantsIterator.hasNext()) {
	            VariantContext variant = variantsIterator.next();
	            printVariantDetails(variant);  // Puedes cambiar este método para procesar la variante como prefieras
	        }
	    }

	    // Método para imprimir los detalles de una variante
	    @SuppressWarnings("deprecation")
		private void printVariantDetails(VariantContext variant) {
	        System.out.print(variant.getChr()+"\t");
	        System.out.print(variant.getStart()+"\t");
	        System.out.print(variant.getReference()+"\t");
	        System.out.print(variant.getAlternateAlleles()+"\t");
	        //System.out.print(variant.getGenotypes()+"\t");
	        
	        for (Genotype genotype : variant.getGenotypes()) {
	        	Map<String, Object> attributes = genotype.getExtendedAttributes();
	        	String BSDP = (String) attributes.get("BSDP");  // Base Depth (BSDP) es un campo común en FORMAT
	            String ACN = (String) attributes.get("ACN"); // Allele Count (ACN)
	            System.out.print("["+genotype.getGenotypeString()+"|"+BSDP+"|"+ACN+"]\t");
	        }
	             
	        System.out.println("");
	    }

	    // Main para probar la clase
	    public static void main(String[] args) {
	        String vcfFilePath = "/Users/estuvar4/Downloads/tmp2.vcf";  // Cambia la ruta al archivo VCF que deseas leer
	        VCF vcfReader = new VCF(vcfFilePath);

	        // Obtener e imprimir el encabezado
	        VCFHeader header = vcfReader.getHeader();
	        System.out.println("Header del VCF:");
	        System.out.println(header);

	        // Iterar sobre las variantes e imprimir detalles
	        System.out.println("\nVariantes:");
	        vcfReader.iterateVariants();
	    }
	}