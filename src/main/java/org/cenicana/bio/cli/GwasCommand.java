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
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "gwas", description = "Perform Polyploid GWAS using EMMAX (P3D) algorithm.")
public class GwasCommand implements Callable<Integer> {

    @Option(names = { "-v", "--vcf" }, description = "Input VCF file", required = true)
    private String vcfPath;

    @Option(names = { "--pheno" }, description = "Phenotype file (CSV/TSV)", required = true)
    private String phenoPath;

    @Option(names = { "--trait" }, description = "Target trait name in phenotype file", required = true)
    private String traitName;

    @Option(names = { "-p", "--ploidy" }, description = "Ploidy level (default: 2)", defaultValue = "2")
    private int ploidy;

    @Option(names = "--loco", description = "Use Leave-One-Chromosome-Out (LOCO) method for better QTL detection")
    private boolean useLoco = false;

    @Option(names = "--epistasis", description = "Run targeted epistasis scan for top hits")
    private boolean runEpistasis = false;

    @Option(names = { "-w",
            "--window" }, description = "Window size for haplotype-block analysis (default: 1 for single SNP)", defaultValue = "1")
    private int windowSize = 1;

    @Option(names = { "-k", "--kinship" }, description = "Path to kinship matrix CSV (Optional, calculated if missing)")
    private String kinshipPath;

    @Option(names = { "-o",
            "--output" }, description = "Output HTML dashboard path", defaultValue = "gwas_results.html")
    private String outputPath;

    @Option(names = { "--qmatrix" }, description = "Path to population structure Q-matrix CSV (Sample, PC1, PC2...)")
    private String qmatrixPath;

    @Option(names = {
            "--ld-prune" }, description = "LD pruning threshold (r2) to reduce marker density (default: 1.0, no pruning)", defaultValue = "1.0")
    private double ldPruneThreshold = 1.0;

    @Option(names = { "--impute" }, description = "Imputation mode: mean, knn (default: mean)", defaultValue = "mean")
    private String imputationMode = "mean";

    @Option(names = { "--maf" }, description = "Min MAF for kinship calculation (default: 0.05)", defaultValue = "0.05")
    private double minMaf;

    @Option(names = { "--fixed" }, description = "Comma-separated list of fixed effects from phenotype file")
    private String fixedEffects;

    @Option(names = {
            "--models" }, description = "Genetic models to test: ADDITIVE, SIMPLEX_DOMINANT, DUPLEX_DOMINANT, TRIPLEX_DOMINANT, SIMPLEX_DOMINANT_REF, DUPLEX_DOMINANT_REF, TRIPLEX_DOMINANT_REF, GENERAL (default: ADDITIVE)", split = ",")
    private List<GwasEngine.GeneticModel> geneticModels = List.of(GwasEngine.GeneticModel.ADDITIVE);

    @Option(names = {
            "--max-missing" }, description = "Maximum allowed missing data per marker (0.0 to 1.0). Default: 0.1 (10%)", defaultValue = "0.1")
    private double maxMissing;

    @Option(names = { "--n-pc" }, description = "Number of principal components to include as covariates (default: 5)", defaultValue = "5")
    private int nPc;

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
            PopulationStructureAnalyzer.PcaResult res = analyzer.computePCA(vcfPath, ploidy, nPc, minMaf, 0.2, true);
            kinship = res.kinshipMatrix;
            pcaCovar = res.pcMatrix;
        }

        // 2.5 Load Q-Matrix if provided
        double[][] qCovar = null;
        if (qmatrixPath != null) {
            System.out.println("[GWAS] Loading Q-Matrix from " + qmatrixPath + "...");
            qCovar = loadQMatrix(qmatrixPath, sampleNames);
        }

        // 3. Run GWAS Engine
        GwasEngine engine = new GwasEngine(ploidy, sampleNames);
        engine.setKinship(kinship);
        engine.setLoco(useLoco);
        engine.setWindowSize(windowSize);
        engine.setRunEpistasis(runEpistasis);
        engine.setModels(new java.util.HashSet<>(geneticModels));
        engine.setLdPruneThreshold(ldPruneThreshold);
        engine.setImputationMode(imputationMode);
        engine.setMaxMissing(maxMissing);
        engine.setMafThreshold(minMaf);

        // Combine PCA + Q-Matrix + Fixed Effects
        double[][] combinedCovar = engine.combineAllCovariates(pcaCovar, qCovar, pheno, fixedCols);
        engine.setCovariates(combinedCovar);

        GwasEngine.GwasResult result = engine.run(vcfPath, yValues, traitName);

        // 4. Generate Dashboard
        System.out.println("[GWAS] Generating interactive dashboard...");
        GwasDashboardGenerator visualizer = new GwasDashboardGenerator();
        visualizer.generate(result, traitName, outputPath, ploidy);

        System.out.println("[GWAS] SUCCESS! Results saved to: " + outputPath);
        return 0;
    }

    private double[][] loadMatrix(String path, int n) throws Exception {
        double[][] matrix = new double[n][n];
        try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(path))) {
            br.readLine(); // Header
            for (int i = 0; i < n; i++) {
                String line = br.readLine();
                if (line == null)
                    break;
                String[] parts = line.split(",");
                for (int j = 0; j < n; j++) {
                    matrix[i][j] = Double.parseDouble(parts[j + 1]);
                }
            }
        }
        return matrix;
    }

    private double[][] loadQMatrix(String path, String[] sampleNames) throws Exception {
        Map<String, double[]> data = new HashMap<>();
        int nCols = 0;
        try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(path))) {
            String header = br.readLine();
            nCols = header.split(",").length - 1;
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(",");
                double[] row = new double[nCols];
                for (int i = 0; i < nCols; i++)
                    row[i] = Double.parseDouble(parts[i + 1]);
                data.put(parts[0], row);
            }
        }
        double[][] qMatrix = new double[sampleNames.length][nCols];
        for (int i = 0; i < sampleNames.length; i++) {
            if (data.containsKey(sampleNames[i])) {
                qMatrix[i] = data.get(sampleNames[i]);
            }
        }
        return qMatrix;
    }
}
