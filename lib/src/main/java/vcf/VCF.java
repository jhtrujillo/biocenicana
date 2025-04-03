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

    /**
     * Constructor that reads the VCF file.
     * @param vcfFilePath Path to the VCF file.
     */
    public VCF(String vcfFilePath) {
        File vcfFile = new File(vcfFilePath);
        VCFFileReader reader = new VCFFileReader(vcfFile, false);

        this.header = reader.getFileHeader();
        this.variantsIterator = reader.iterator();
    }

    /**
     * Method to get the header.
     * @return The VCF file header.
     */
    public VCFHeader getHeader() {
        return header;
    }

    /**
     * Method to set the header.
     * @param header The VCF file header.
     */
    public void setHeader(VCFHeader header) {
        this.header = header;
    }

    /**
     * Method to get the next VariantContext (a variant).
     * @return The next VariantContext, or null if there are no more variants.
     */
    public VariantContext getNextVariant() {
        if (variantsIterator.hasNext()) {
            return variantsIterator.next();
        }
        return null;  // If there are no more variants
    }

    /**
     * Method to iterate over all variants.
     */
    public void iterateVariants() {
        while (variantsIterator.hasNext()) {
            VariantContext variant = variantsIterator.next();
            printVariantDetails(variant);  // You can change this method to process the variant as you prefer
        }
    }

    /**
     * Method to print the details of a variant.
     * @param variant The VariantContext to print.
     */
    @SuppressWarnings("deprecation")
    private void printVariantDetails(VariantContext variant) {
        System.out.print(variant.getChr() + "\t");
        System.out.print(variant.getStart() + "\t");
        System.out.print(variant.getReference() + "\t");
        System.out.print(variant.getAlternateAlleles() + "\t");

        for (Genotype genotype : variant.getGenotypes()) {
            Map<String, Object> attributes = genotype.getExtendedAttributes();
            String BSDP = (String) attributes.get("BSDP");  // Base Depth (BSDP) is a common field in FORMAT
            String ACN = (String) attributes.get("ACN"); // Allele Count (ACN)
            System.out.print("[" + genotype.getGenotypeString() + "|" + BSDP + "|" + ACN + "]\t");
        }

        System.out.println("");
    }

    /**
     * Main method to test the class.
     * @param args Command line arguments.
     */
    public static void main(String[] args) {
        String vcfFilePath = "/Users/estuvar4/Downloads/tmp2.vcf";  // Change the path to the VCF file you want to read
        VCF vcfReader = new VCF(vcfFilePath);

        // Get and print the header
        VCFHeader header = vcfReader.getHeader();
        System.out.println("VCF Header:");
        System.out.println(header);

        // Iterate over the variants and print details
        System.out.println("\nVariants:");
        vcfReader.iterateVariants();
    }
}
