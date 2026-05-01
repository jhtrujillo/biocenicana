package org.cenicana.bio.cli;

import org.cenicana.bio.core.PopulationStructureAnalyzer;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.io.File;
import java.util.concurrent.Callable;

@Command(name = "pop-structure", description = "Calculate population structure using PCA (Principal Component Analysis)")
public class PopStructureCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, required = true, description = "Input VCF file")
    private String vcfFile;

    @Option(names = {"-o", "--output"}, required = true, description = "Output base name (e.g., 'results' creates results.pca.csv)")
    private String outputBase;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level (e.g., 10 for sugarcane)", defaultValue = "10")
    private int ploidy;

    @Option(names = {"--pcs"}, description = "Number of Principal Components to compute", defaultValue = "10")
    private int numPCs;

    @Option(names = {"--min-maf"}, description = "Minimum MAF for SNP filtering", defaultValue = "0.05")
    private double minMaf;

    @Option(names = {"--max-missing"}, description = "Maximum missingness per SNP", defaultValue = "0.2")
    private double maxMissing;

    @Override
    public Integer call() throws Exception {
        File f = new File(vcfFile);
        if (!f.exists()) {
            System.err.println("Error: VCF file not found: " + vcfFile);
            return 1;
        }

        System.out.println("=================================================");
        System.out.println("BioJava Population Structure (PCA)");
        System.out.println("=================================================");
        System.out.println("Input:      " + vcfFile);
        System.out.println("Output:     " + outputBase + ".pca.csv");
        System.out.println("Ploidy:     " + ploidy);
        System.out.println("PCs:        " + numPCs);
        System.out.println("Filters:    MAF > " + minMaf + ", Missing < " + maxMissing);
        System.out.println("=================================================\n");

        PopulationStructureAnalyzer analyzer = new PopulationStructureAnalyzer();
        PopulationStructureAnalyzer.PcaResult result = analyzer.computePCA(vcfFile, ploidy, numPCs, minMaf, maxMissing);
        
        analyzer.exportPca(result, outputBase + ".pca.csv");
        analyzer.exportKinship(result, outputBase + ".kinship.csv");
        
        String htmlPath = outputBase + ".pca.html";
        org.cenicana.bio.io.PcaDashboardGenerator.generateReport(result, htmlPath);

        System.out.println("\n✅  PCA analysis complete.");
        
        // Report explained variance
        System.out.println("\nExplained Variance per PC:");
        for (int i = 0; i < result.explainedVariance.length; i++) {
            System.out.printf(" - PC%d: %.2f%%%n", (i+1), result.explainedVariance[i] * 100.0);
        }

        return 0;
    }
}
