package org.cenicana.bio.cli;

import org.cenicana.bio.core.SnpClusterAnalyzer;
import org.cenicana.bio.core.SnpClusterAnalyzer.SnpResult;
import org.cenicana.bio.io.SnpExplorerDashboard;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.List;
import java.util.concurrent.Callable;

@Command(name = "snp-explorer", 
         mixinStandardHelpOptions = true, 
         description = "Generates an interactive HTML dashboard to explore dosage groupings for each SNP.")
public class SnpExplorerCommand implements Callable<Integer> {

    @Option(names = { "-v", "--vcf" }, description = "Path to the input VCF file (required for AD-Plots).")
    private String vcfFile;

    @Option(names = { "-m", "--matrix" }, description = "Path to the dosage matrix TSV (optional if VCF is provided).")
    private String matrixFile;

    @Option(names = { "-p", "--ploidy" }, required = true, description = "Ploidy level of the population.")
    private int ploidy;

    @Option(names = { "--pca" }, description = "Optionally, path to the PCA results CSV (from pop-structure).")
    private String pcaFile;

    @Option(names = { "-o", "--output" }, defaultValue = "snp_explorer.html", description = "Output HTML file path.")
    private String outputFile;

    @Override
    public Integer call() throws Exception {
        if (vcfFile == null && matrixFile == null) {
            System.err.println("Error: You must provide either --vcf or --matrix.");
            return 1;
        }

        SnpClusterAnalyzer analyzer = new SnpClusterAnalyzer();
        List<SnpClusterAnalyzer.SnpResult> results;

        if (vcfFile != null) {
            System.out.println("[SNP-Explorer] Reading VCF (extracting dosages + AD): " + vcfFile);
            results = analyzer.analyzeVcf(vcfFile, ploidy);
        } else {
            System.out.println("[SNP-Explorer] Reading matrix: " + matrixFile);
            results = analyzer.analyzeMatrix(matrixFile, ploidy);
        }
        
        List<SnpClusterAnalyzer.SampleCoord> coords = null;
        if (pcaFile != null) {
            System.out.println("[SNP-Explorer] Loading PCA coordinates: " + pcaFile);
            coords = analyzer.loadPcaCsv(pcaFile);
        }
        
        System.out.println("[SNP-Explorer] Processed " + results.size() + " SNPs.");
        System.out.println("[SNP-Explorer] Generating dashboard: " + outputFile);
        
        SnpExplorerDashboard.generate(results, coords, ploidy, outputFile);
        
        System.out.println("\n✅  SNP Explorer dashboard ready: " + outputFile);
        return 0;
    }
}
