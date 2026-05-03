package org.cenicana.bio.cli;

import org.cenicana.bio.core.GwasEngine;
import org.cenicana.bio.core.PopulationStructureAnalyzer;
import org.cenicana.bio.io.PhenotypeData;
import org.cenicana.bio.io.VcfFastReader;
import org.cenicana.bio.io.GwasDashboardGenerator;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

@Command(name = "gwas", description = "Perform Polyploid GWAS using EMMAX (P3D) algorithm.")
public class GwasCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, description = "Input VCF file", required = true)
    private String vcfPath;

    @Option(names = {"--pheno"}, description = "Phenotype file (CSV/TSV)", required = true)
    private String phenoPath;

    @Option(names = {"--trait"}, description = "Target trait name in phenotype file", required = true)
    private String traitName;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level (default: 2)", defaultValue = "2")
    private int ploidy;

    @Option(names = "--loco", description = "Use Leave-One-Chromosome-Out (LOCO) method for better QTL detection")
    private boolean useLoco = false;

    @Option(names = {"-k", "--kinship"}, description = "Path to kinship matrix CSV (Optional, calculated if missing)")
    private String kinshipPath;

    @Option(names = {"-o", "--output"}, description = "Output HTML dashboard path", defaultValue = "gwas_results.html")
    private String outputPath;

    @Option(names = {"--maf"}, description = "Min MAF for kinship calculation (default: 0.05)", defaultValue = "0.05")
    private double minMaf;

    @Option(names = {"--fixed"}, description = "Comma-separated list of fixed effects from phenotype file")
    private String fixedEffects;

    @Override
    public Integer call() throws Exception {
        System.out.println("\n[BioJava] Starting Polyploid GWAS Pipeline...");
        
        // 1. Load Phenotypes
        List<String> fixedCols = fixedEffects != null ? List.of(fixedEffects.split(",")) : null;
        PhenotypeData pheno = new PhenotypeData();
        pheno.load(phenoPath, traitName, fixedCols);
        
        String[] sampleNames = VcfFastReader.getSampleIds(vcfPath);
        double[] yValues = pheno.getOrderedValues(sampleNames);
        
        // 2. Load or Compute Kinship and PCA
        double[][] kinship;
        double[][] pcaCovar = null;
        if (kinshipPath != null && new File(kinshipPath).exists()) {
            System.out.println("[GWAS] Loading Kinship from " + kinshipPath + "...");
            kinship = loadMatrix(kinshipPath, sampleNames.length);
        } else {
            System.out.println("[GWAS] Kinship not provided. Computing VanRaden matrix from VCF...");
            PopulationStructureAnalyzer analyzer = new PopulationStructureAnalyzer();
            PopulationStructureAnalyzer.PcaResult res = analyzer.computePCA(vcfPath, ploidy, 5, minMaf, 0.2, true);
            kinship = res.kinshipMatrix;
            pcaCovar = res.pcMatrix;
        }

        // 3. Run GWAS Engine
        GwasEngine engine = new GwasEngine(ploidy, sampleNames);
        engine.setKinship(kinship);
        engine.setLoco(useLoco);
        
        // Combine PCA + Fixed Effects
        if (pcaCovar != null || fixedCols != null) {
            double[][] combinedCovar = engine.combineCovariates(pcaCovar, pheno, fixedCols);
            engine.setCovariates(combinedCovar);
        }
        
        List<GwasEngine.GwasHit> hits = engine.run(vcfPath, yValues, traitName);

        // 4. Generate Dashboard
        System.out.println("[GWAS] Generating interactive dashboard...");
        GwasDashboardGenerator visualizer = new GwasDashboardGenerator();
        visualizer.generate(hits, traitName, outputPath, ploidy);

        System.out.println("[GWAS] SUCCESS! Results saved to: " + outputPath);
        return 0;
    }

    private double[][] loadMatrix(String path, int n) throws Exception {
        double[][] matrix = new double[n][n];
        try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(path))) {
            br.readLine(); // Header
            for (int i = 0; i < n; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] parts = line.split(",");
                for (int j = 0; j < n; j++) {
                    matrix[i][j] = Double.parseDouble(parts[j + 1]);
                }
            }
        }
        return matrix;
    }
}
