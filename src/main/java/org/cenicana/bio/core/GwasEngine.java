package org.cenicana.bio.core;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import org.cenicana.bio.utils.GwasMathUtils;
import org.cenicana.bio.io.PhenotypeData;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Modernized GWAS Engine supporting Mixed Models (EMMAX/P3D), 
 * Polyploidy, Fixed Effects, Partial R2, and LOCO (Leave-One-Chromosome-Out).
 */
public class GwasEngine {

    public enum GeneticModel {
        ADDITIVE, SIMPLEX_DOMINANT, DUPLEX_DOMINANT, 
        SIMPLEX_DOMINANT_REF, DUPLEX_DOMINANT_REF, GENERAL
    }

    public static class GwasHit {
        public String markerId;
        public String chromosome;
        public long position;
        public String refAllele;
        public String altAllele;
        public double pValue;
        public double qValue; // FDR
        public double effect;
        public double r2;
        public double aic;
        public String model; // Additive, SimplexDominant, etc.
        public boolean isBestModel; // To flag the best model for this SNP
        public List<Double>[] phenotypesByDosage; // To store distributions for boxplots
        public List<String>[] samplesByDosage; // To track elite candidates
    }

    public static class GwasResult {
        public List<GwasHit> hits;
        public List<GwasInteraction> interactions = new ArrayList<>();
    }

    public static class GwasInteraction {
        public String marker1;
        public String marker2;
        public double pValue;
        public double effect;
    }

    private static class MarkerData {
        String id;
        String chrom;
        long pos;
        String ref;
        String alt;
        double[] dosages;
    }

    private int ploidy;
    private String[] sampleNames;
    private double[][] covariates; // Population structure (PCA) + Fixed effects
    private double[][] kinship;
    private boolean useLoco = false;
    private boolean runEpistasis = false;
    private int windowSize = 1;
    private int numThreads = Runtime.getRuntime().availableProcessors();
    private Set<GeneticModel> models = new HashSet<>(List.of(GeneticModel.ADDITIVE));
    private double ldPruneThreshold = 1.0;
    private double maxMissing = 0.1;
    private double mafThreshold = 0.01;
    private String imputationMode = "mean"; // mean or knn

    public GwasEngine(int ploidy, String[] sampleNames) {
        this.ploidy = ploidy;
        this.sampleNames = sampleNames;
    }

    public void setKinship(double[][] kinship) {
        this.kinship = kinship;
    }

    public void setLoco(boolean useLoco) {
        this.useLoco = useLoco;
    }

    public void setRunEpistasis(boolean runEpistasis) {
        this.runEpistasis = runEpistasis;
    }

    public void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
    }

    public void setModels(Set<GeneticModel> models) {
        this.models = models;
    }

    public void setLdPruneThreshold(double ldPruneThreshold) {
        this.ldPruneThreshold = ldPruneThreshold;
    }

    public void setImputationMode(String mode) {
        this.imputationMode = mode;
    }

    public void setMaxMissing(double maxMissing) {
        this.maxMissing = maxMissing;
    }

    public void setMafThreshold(double mafThreshold) {
        this.mafThreshold = mafThreshold;
    }

    public void setCovariates(double[][] covariates) { this.covariates = covariates; }
    public void setNumThreads(int numThreads) { this.numThreads = numThreads; }

    public double[][] combineAllCovariates(double[][] pcaCovar, double[][] qCovar, PhenotypeData pheno, List<String> fixedCols) {
        int n = sampleNames.length;
        int pcaCount = (pcaCovar != null) ? pcaCovar[0].length : 0;
        int qCount = (qCovar != null) ? qCovar[0].length : 0;
        
        List<double[]> dummyEncoded = new ArrayList<>();
        if (fixedCols != null) {
            for (String colName : fixedCols) {
                Map<String, List<Integer>> levelIndices = new HashMap<>();
                for (int i = 0; i < n; i++) {
                    String val = pheno.getFixedEffects(sampleNames[i]).get(colName);
                    if (val != null) levelIndices.computeIfAbsent(val, k -> new ArrayList<>()).add(i);
                }
                if (levelIndices.size() < 2) continue;
                List<String> uniqueLevels = new ArrayList<>(levelIndices.keySet());
                for (int i = 0; i < uniqueLevels.size() - 1; i++) {
                    double[] dummy = new double[n];
                    for (int idx : levelIndices.get(uniqueLevels.get(i))) dummy[idx] = 1.0;
                    dummyEncoded.add(dummy);
                }
            }
        }
        
        int totalCovs = pcaCount + qCount + dummyEncoded.size();
        if (totalCovs == 0) return null;
        
        double[][] combined = new double[n][totalCovs];
        for (int i = 0; i < n; i++) {
            int currentPos = 0;
            if (pcaCovar != null) {
                System.arraycopy(pcaCovar[i], 0, combined[i], currentPos, pcaCount);
                currentPos += pcaCount;
            }
            if (qCovar != null) {
                System.arraycopy(qCovar[i], 0, combined[i], currentPos, qCount);
                currentPos += qCount;
            }
            for (int j = 0; j < dummyEncoded.size(); j++) {
                combined[i][currentPos + j] = dummyEncoded.get(j)[i];
            }
        }
        return combined;
    }

    public GwasResult run(String vcfPath, double[] yValues, String traitName) throws Exception {
        System.out.println("[GWAS] Loading variants and filtering...");
        
        // 1. Filter individuals with missing phenotypes
        List<Integer> validIndices = new ArrayList<>();
        for (int i = 0; i < yValues.length; i++) {
            if (!Double.isNaN(yValues[i])) validIndices.add(i);
        }
        int nFiltered = validIndices.size();
        double[] Yf = new double[nFiltered];
        for (int i = 0; i < nFiltered; i++) Yf[i] = yValues[validIndices.get(i)];
        
        // 2. Load and Filter Markers into memory
        Map<String, List<MarkerData>> chromToMarkers = new LinkedHashMap<>();
        int totalLoaded = 0;
        
        try (BufferedReader br = new BufferedReader(new FileReader(vcfPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] cols = line.split("\t");
                double[] dosages = extractDosages(cols, validIndices, kinship);
                if (dosages == null) continue;

                MarkerData md = new MarkerData();
                md.id = cols[2].equals(".") ? cols[0] + ":" + cols[1] : cols[2];
                md.chrom = cols[0];
                md.pos = Long.parseLong(cols[1]);
                md.ref = cols[3];
                md.alt = cols[4];
                md.dosages = dosages;

                chromToMarkers.computeIfAbsent(md.chrom, k -> new ArrayList<>()).add(md);
                totalLoaded++;
            }
        }
        
        // 2.2 LD Pruning if requested
        if (ldPruneThreshold < 0.999) {
            System.out.println("[GWAS] Pruning markers by LD (r2 < " + ldPruneThreshold + ")...");
            int before = totalLoaded;
            pruneByLD(chromToMarkers);
            totalLoaded = chromToMarkers.values().stream().mapToInt(List::size).sum();
            System.out.println("[GWAS] Pruning complete: " + before + " -> " + totalLoaded + " markers.");
        }

        // Smart Grouping for LOCO: Group chroms with < 100 markers into a single "pseudo-chrom"
        Map<String, List<MarkerData>> locoGroups = new LinkedHashMap<>();
        if (useLoco) {
            List<MarkerData> smallContigs = new ArrayList<>();
            for (Map.Entry<String, List<MarkerData>> entry : chromToMarkers.entrySet()) {
                if (entry.getValue().size() < 100) {
                    smallContigs.addAll(entry.getValue());
                } else {
                    locoGroups.put(entry.getKey(), entry.getValue());
                }
            }
            if (!smallContigs.isEmpty()) {
                locoGroups.put("__SMALL_CONTIGS__", smallContigs);
            }
        } else {
            locoGroups = chromToMarkers;
        }

        System.out.println("[GWAS] Total markers: " + totalLoaded + " in " + chromToMarkers.size() + " sequences.");
        if (useLoco) System.out.println("[GWAS] LOCO Groups formed: " + locoGroups.size());

        // 3. Prepare Covariates (W)
        int numCovs = (covariates != null ? covariates[0].length : 0);
        DMatrixRMaj W = new DMatrixRMaj(nFiltered, 1 + numCovs);
        for (int i = 0; i < nFiltered; i++) W.set(i, 0, 1.0);
        if (covariates != null) {
            for (int i = 0; i < nFiltered; i++) {
                int originalIdx = validIndices.get(i);
                for (int j = 0; j < numCovs; j++) {
                    W.set(i, j + 1, covariates[originalIdx][j]);
                }
            }
        }

        List<GwasHit> allHits = Collections.synchronizedList(new ArrayList<>());
        String[] filteredNames = new String[nFiltered];
        for (int i = 0; i < nFiltered; i++) filteredNames[i] = sampleNames[validIndices.get(i)];

        if (!useLoco) {
            // GLOBAL MODE
            System.out.println("[GWAS] Running Global Analysis...");
            if (kinship == null) {
                List<double[]> allDosages = new ArrayList<>();
                for (List<MarkerData> list : chromToMarkers.values()) {
                    for (MarkerData md : list) allDosages.add(md.dosages);
                }
                kinship = PopulationStructureAnalyzer.calculateVanRadenKinship(allDosages, ploidy);
            }
            runChromBlock(chromToMarkers.keySet(), chromToMarkers, kinship, nFiltered, filteredNames, Yf, W, allHits);
        } else {
            // LOCO MODE
            System.out.println("[GWAS] Running LOCO Analysis (Leave-One-Chromosome-Out)...");
            for (String groupName : locoGroups.keySet()) {
                System.out.println("[GWAS] Analyzing Group: " + groupName + " (LOCO)");
                
                // Build Kinship excluding markers in THIS group
                List<double[]> dosagesExcl = new ArrayList<>();
                for (Map.Entry<String, List<MarkerData>> entry : locoGroups.entrySet()) {
                    if (!entry.getKey().equals(groupName)) {
                        for (MarkerData md : entry.getValue()) dosagesExcl.add(md.dosages);
                    }
                }
                
                double[][] locoKinship = PopulationStructureAnalyzer.calculateVanRadenKinship(dosagesExcl, ploidy);
                
                // Scan the markers in this group (could be one chrom or many small contigs)
                Map<String, List<MarkerData>> markersInGroup = new HashMap<>();
                if (groupName.equals("__SMALL_CONTIGS__")) {
                    for (MarkerData md : locoGroups.get(groupName)) {
                        markersInGroup.computeIfAbsent(md.chrom, k -> new ArrayList<>()).add(md);
                    }
                } else {
                    markersInGroup.put(groupName, locoGroups.get(groupName));
                }
                
                runChromBlock(markersInGroup.keySet(), markersInGroup, locoKinship, nFiltered, filteredNames, Yf, W, allHits);
            }
        }

        allHits.sort((a, b) -> Double.compare(a.pValue, b.pValue));
        
        // 4. Apply FDR (Benjamini-Hochberg)
        applyFDRCorrection(allHits);
        
        // 5. Select Best Model per Marker (based on AIC)
        selectBestModels(allHits);

        GwasResult result = new GwasResult();
        result.hits = allHits;

        if (runEpistasis) {
            // We need weights, U, Ystar, Wstar from the GLOBAL model for epistasis leads
            // For simplicity in LOCO, we re-run a global kinship/decomposition if needed or use the last one
            // Let's use the parameters from a global run
            System.out.println("[GWAS] Running epistasis phase...");
            List<double[]> allDosages = new ArrayList<>();
            for (List<MarkerData> list : chromToMarkers.values()) {
                for (MarkerData md : list) allDosages.add(md.dosages);
            }
            double[][] globalK = PopulationStructureAnalyzer.calculateVanRadenKinship(allDosages, ploidy);
            
            // Temporary block to get global decomposition
            DMatrixRMaj K = new DMatrixRMaj(globalK);
            EigenDecomposition_F64<DMatrixRMaj> evd = DecompositionFactory_DDRM.eig(nFiltered, true);
            evd.decompose(K);
            double[] lambdas = new double[nFiltered];
            DMatrixRMaj U = new DMatrixRMaj(nFiltered, nFiltered);
            for (int i = 0; i < nFiltered; i++) {
                lambdas[i] = evd.getEigenvalue(i).getReal();
                DMatrixRMaj v = evd.getEigenVector(i);
                for (int j = 0; j < nFiltered; j++) U.set(j, i, v.get(j, 0));
            }
            DMatrixRMaj Ymat = new DMatrixRMaj(nFiltered, 1);
            for (int i = 0; i < nFiltered; i++) Ymat.set(i, 0, Yf[i]);
            DMatrixRMaj Ystar = new DMatrixRMaj(nFiltered, 1);
            CommonOps_DDRM.multTransA(U, Ymat, Ystar);
            DMatrixRMaj Wstar = new DMatrixRMaj(nFiltered, W.numCols);
            CommonOps_DDRM.multTransA(U, W, Wstar);
            double delta = 1.0; 
            double[] weights = new double[nFiltered];
            for (int i = 0; i < nFiltered; i++) weights[i] = 1.0 / (Math.max(1e-6, lambdas[i]) + delta);

            result.interactions = runEpistasisScan(allHits, chromToMarkers, Yf, W, weights, U, Ystar, Wstar);
        }

        return result;
    }

    private void runChromBlock(Set<String> chromsToScan, Map<String, List<MarkerData>> chromToMarkers, double[][] kMatrix, 
                               int nFiltered, String[] filteredNames, double[] Yf, DMatrixRMaj W, List<GwasHit> allHits) {
        // 1. EVD of Kinship
        DMatrixRMaj K = new DMatrixRMaj(kMatrix);
        EigenDecomposition_F64<DMatrixRMaj> evd = DecompositionFactory_DDRM.eig(nFiltered, true);
        if (!evd.decompose(K)) {
            System.err.println("[GWAS] ERROR: Kinship decomposition failed.");
            return;
        }

        double[] lambdas = new double[nFiltered];
        DMatrixRMaj U = new DMatrixRMaj(nFiltered, nFiltered);
        for (int i = 0; i < nFiltered; i++) {
            lambdas[i] = evd.getEigenvalue(i).getReal();
            DMatrixRMaj v = evd.getEigenVector(i);
            for (int j = 0; j < nFiltered; j++) U.set(j, i, v.get(j, 0));
        }

        // 2. Rotate Y and W: Y* = U'Y, W* = U'W
        DMatrixRMaj Ymat = new DMatrixRMaj(nFiltered, 1);
        for (int i = 0; i < nFiltered; i++) Ymat.set(i, 0, Yf[i]);
        DMatrixRMaj Ystar = new DMatrixRMaj(nFiltered, 1);
        CommonOps_DDRM.multTransA(U, Ymat, Ystar);

        DMatrixRMaj Wstar = new DMatrixRMaj(nFiltered, W.numCols);
        CommonOps_DDRM.multTransA(U, W, Wstar);

        double delta = 1.0; 
        double[] weights = new double[nFiltered];
        for (int i = 0; i < nFiltered; i++) weights[i] = 1.0 / (Math.max(1e-6, lambdas[i]) + delta);
        final double rssReduced = calculateRSS(Ystar, Wstar, lambdas, delta);

        // 3. Parallel Scan
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        for (String chrom : chromsToScan) {
            List<MarkerData> markers = chromToMarkers.get(chrom);
            if (markers == null) continue;
            
            // Sliding Window grouping
            for (int i = 0; i < markers.size(); i += windowSize) {
                final List<MarkerData> window = markers.subList(i, Math.min(i + windowSize, markers.size()));
                executor.submit(() -> {
                    int q = Wstar.numCols;
                    DMatrixRMaj M = new DMatrixRMaj(nFiltered, q + 1);
                    DMatrixRMaj MtVinv = new DMatrixRMaj(q + 1, nFiltered);
                    DMatrixRMaj LHS = new DMatrixRMaj(q + 1, q + 1);
                    DMatrixRMaj RHS = new DMatrixRMaj(q + 1, 1);
                    DMatrixRMaj beta = new DMatrixRMaj(q + 1, 1);
                    DMatrixRMaj invLHS = new DMatrixRMaj(q + 1, q + 1);
                    DMatrixRMaj Xmat = new DMatrixRMaj(nFiltered, 1);
                    DMatrixRMaj Xstar = new DMatrixRMaj(nFiltered, 1);

                    if (window.size() == 1) {
                        MarkerData md = window.get(0);
                        for (GeneticModel model : models) {
                            double[] X = recodeForModel(md.dosages, model);
                            if (X == null) continue;
                            
                            int df = (model == GeneticModel.GENERAL) ? getUniqueDosageCount(md.dosages) - 1 : 1;
                            GwasHit hit = testModel(md.id, md.chrom, md.pos, md.ref, md.alt, X, df, filteredNames, Yf, Ystar, Wstar, U, weights, rssReduced, model.toString(), M, MtVinv, LHS, RHS, beta, invLHS, Xmat, Xstar);
                            if (hit.pValue < 1.0) allHits.add(hit);
                        }
                    } else {
                        double[] haploDosage = computeHaploDosage(window);
                        MarkerData first = window.get(0);
                        MarkerData last = window.get(window.size() - 1);
                        String blockId = String.format("Block_%s:%d-%d", first.chrom, first.pos, last.pos);
                        GwasHit hit = testModel(blockId, first.chrom, (first.pos + last.pos)/2, "Block", "Block", haploDosage, 1, filteredNames, Yf, Ystar, Wstar, U, weights, rssReduced, "Haplotype", M, MtVinv, LHS, RHS, beta, invLHS, Xmat, Xstar);
                        hit.markerId += " (" + window.size() + " SNPs)";
                        allHits.add(hit);
                    }
                });
            }
        }
        executor.shutdown();
        try { executor.awaitTermination(1, TimeUnit.HOURS); } catch (InterruptedException e) { e.printStackTrace(); }
    }

    private double[] computeHaploDosage(List<MarkerData> window) {
        int n = window.get(0).dosages.length;
        int m = window.size();
        DMatrixRMaj mat = new DMatrixRMaj(n, m);
        for (int j = 0; j < m; j++) {
            double[] d = window.get(j).dosages;
            for (int i = 0; i < n; i++) mat.set(i, j, d[i]);
        }
        
        // Center the matrix for SVD
        for (int j = 0; j < m; j++) {
            double sum = 0;
            for (int i = 0; i < n; i++) sum += mat.get(i, j);
            double mean = sum / n;
            for (int i = 0; i < n; i++) mat.set(i, j, mat.get(i, j) - mean);
        }
        
        SingularValueDecomposition_F64<DMatrixRMaj> svd = DecompositionFactory_DDRM.svd(n, m, true, true, false);
        if (!svd.decompose(mat)) return window.get(0).dosages;
        
        DMatrixRMaj U = svd.getU(null, false);
        double[] pc1 = new double[n];
        for (int i = 0; i < n; i++) pc1[i] = U.get(i, 0);
        
        // Scale to 0..ploidy range to maintain consistency with boxplots
        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
        for (double v : pc1) { min = Math.min(min, v); max = Math.max(max, v); }
        for (int i = 0; i < n; i++) pc1[i] = ((pc1[i] - min) / (max - min + 1e-10)) * ploidy;
        
        return pc1;
    }

    private double calculateRSS(DMatrixRMaj Ystar, DMatrixRMaj Wstar, double[] lambdas, double d) {
        int n = Ystar.numRows;
        int q = Wstar.numCols;
        DMatrixRMaj LHS = new DMatrixRMaj(q, q);
        DMatrixRMaj RHS = new DMatrixRMaj(q, 1);
        for (int i = 0; i < n; i++) {
            double w = 1.0 / (lambdas[i] + d);
            for (int j = 0; j < q; j++) {
                RHS.add(j, 0, Wstar.get(i, j) * w * Ystar.get(i, 0));
                for (int k = 0; k < q; k++) LHS.add(j, k, Wstar.get(i, j) * w * Wstar.get(i, k));
            }
        }
        DMatrixRMaj beta = new DMatrixRMaj(q, 1);
        if (!CommonOps_DDRM.solve(LHS, RHS, beta)) return 1e10;
        double rss = 0;
        for (int i = 0; i < n; i++) {
            double pred = 0;
            for (int j = 0; j < q; j++) pred += Wstar.get(i, j) * beta.get(j, 0);
            double err = Ystar.get(i, 0) - pred;
            rss += (err * err) / (lambdas[i] + d);
        }
        return rss;
    }

    private GwasHit testModel(String id, String chrom, long pos, String ref, String alt, double[] X, int df, String[] filteredNames, double[] Yf, DMatrixRMaj Ystar, DMatrixRMaj Wstar, DMatrixRMaj U, double[] weights, double rssReduced, String modelName,
                              DMatrixRMaj M, DMatrixRMaj MtVinv, DMatrixRMaj LHS, DMatrixRMaj RHS, DMatrixRMaj beta, DMatrixRMaj invLHS, DMatrixRMaj Xmat, DMatrixRMaj Xstar) {
        int n = Yf.length;
        int q = Wstar.numCols;
        
        // Handle multi-df (General model)
        int markerCols = df;
        DMatrixRMaj M_multi = new DMatrixRMaj(n, q + markerCols);
        CommonOps_DDRM.insert(Wstar, M_multi, 0, 0);
        
        for (int d = 0; d < markerCols; d++) {
            for (int i = 0; i < n; i++) Xmat.set(i, 0, X[i * markerCols + d]);
            CommonOps_DDRM.multTransA(U, Xmat, Xstar);
            CommonOps_DDRM.insert(Xstar, M_multi, 0, q + d);
        }

        DMatrixRMaj LHS_multi = new DMatrixRMaj(q + markerCols, q + markerCols);
        DMatrixRMaj RHS_multi = new DMatrixRMaj(q + markerCols, 1);
        DMatrixRMaj beta_multi = new DMatrixRMaj(q + markerCols, 1);

        for (int i = 0; i < n; i++) {
            double w = weights[i];
            for (int j = 0; j < q + markerCols; j++) {
                double val = M_multi.get(i, j) * w;
                RHS_multi.add(j, 0, val * Ystar.get(i, 0));
                for (int k = 0; k < q + markerCols; k++) LHS_multi.add(j, k, val * M_multi.get(i, k));
            }
        }

        if (!CommonOps_DDRM.solve(LHS_multi, RHS_multi, beta_multi)) {
            GwasHit h = new GwasHit(); h.pValue = 1.0; return h;
        }

        double rssFull = 0;
        for (int i = 0; i < n; i++) {
            double pred = 0;
            for (int j = 0; j < q + markerCols; j++) pred += M_multi.get(i, j) * beta_multi.get(j, 0);
            rssFull += Math.pow(Ystar.get(i, 0) - pred, 2) * weights[i];
        }

        double partialR2 = (rssReduced - rssFull) / (rssReduced + 1e-10);
        
        // Likelihood Ratio Test or F-test
        double fStat = ((rssReduced - rssFull) / df) / (rssFull / (n - q - df));
        double p = 1.0 - GwasMathUtils.fCDF(fStat, df, n - q - df);
        
        // AIC Calculation: n * ln(RSS/n) + 2k
        double aic = n * Math.log(rssFull / n) + 2 * (q + df + 1);

        GwasHit hit = new GwasHit();
        hit.markerId = id; hit.chromosome = chrom; hit.position = pos;
        hit.refAllele = ref; hit.altAllele = alt;
        hit.effect = beta_multi.get(q, 0); // Primary effect
        hit.model = modelName; hit.r2 = Math.max(0, partialR2);
        hit.pValue = Math.max(p, 1e-300);
        hit.aic = aic;
        
        if (hit.pValue < 0.001) {
            hit.phenotypesByDosage = new List[ploidy + 1];
            hit.samplesByDosage = new List[ploidy + 1];
            for (int i = 0; i <= ploidy; i++) {
                hit.phenotypesByDosage[i] = new ArrayList<>();
                hit.samplesByDosage[i] = new ArrayList<>();
            }
            for (int i = 0; i < n; i++) {
                int dosage = (int) Math.round(X[i * markerCols]); // Use first col for dosage tracking
                if (dosage >= 0 && dosage <= ploidy) {
                    hit.phenotypesByDosage[dosage].add(Yf[i]);
                    hit.samplesByDosage[dosage].add(filteredNames[i]);
                }
            }
        }
        return hit;
    }

    private double[] recodeForModel(double[] dosages, GeneticModel model) {
        int n = dosages.length;
        switch (model) {
            case ADDITIVE: return dosages;
            case SIMPLEX_DOMINANT: return recode(dosages, 1.0);
            case DUPLEX_DOMINANT: return recode(dosages, 2.0);
            case SIMPLEX_DOMINANT_REF: return recodeRef(dosages, ploidy - 1);
            case DUPLEX_DOMINANT_REF: return recodeRef(dosages, ploidy - 2);
            case GENERAL:
                Set<Integer> unique = new TreeSet<>();
                for (double d : dosages) unique.add((int)Math.round(d));
                if (unique.size() < 2) return null;
                List<Integer> levels = new ArrayList<>(unique);
                int df = levels.size() - 1;
                double[] Xmulti = new double[n * df];
                for (int i = 0; i < n; i++) {
                    int d = (int)Math.round(dosages[i]);
                    for (int j = 0; j < df; j++) {
                        Xmulti[i * df + j] = (d == levels.get(j)) ? 1.0 : 0.0;
                    }
                }
                return Xmulti;
            default: return dosages;
        }
    }

    private int getUniqueDosageCount(double[] dosages) {
        Set<Integer> unique = new HashSet<>();
        for (double d : dosages) unique.add((int)Math.round(d));
        return unique.size();
    }

    private double[] extractDosages(String[] cols, List<Integer> validIndices, double[][] kinshipMatrix) {
        int n = sampleNames.length;
        int nFiltered = validIndices.size();
        double[] dosages = new double[nFiltered];
        List<Integer> missingIndices = new ArrayList<>();
        double sum = 0;
        int count = 0;

        for (int i = 0; i < nFiltered; i++) {
            int idx = validIndices.get(i);
            String gData = cols[9 + idx];
            if (gData.startsWith(".")) {
                missingIndices.add(i);
                dosages[i] = Double.NaN;
            } else {
                String gt = gData.split(":")[0];
                int dosage = 0;
                for (char c : gt.toCharArray()) if (c > '0' && c <= '9') dosage++;
                dosages[i] = dosage;
                sum += dosage;
                count++;
            }
        }

        double missingRatio = (double)missingIndices.size() / nFiltered;
        if (count == 0 || missingRatio > maxMissing) return null;
        if (maxMissing == 0 && !missingIndices.isEmpty()) return null;
        double mean = sum / count;
        
        double freq = sum / (count * ploidy);
        double mafValue = Math.min(freq, 1.0 - freq);
        if (mafValue < mafThreshold) return null;

        // Imputation
        for (int i : missingIndices) {
            if ("knn".equalsIgnoreCase(imputationMode) && kinshipMatrix != null) {
                dosages[i] = imputeWithKNN(i, dosages, validIndices, kinshipMatrix, 5);
            } else {
                dosages[i] = mean;
            }
        }
        return dosages;
    }

    private double imputeWithKNN(int targetIdx, double[] dosages, List<Integer> validIndices, double[][] kinship, int k) {
        // Find top k relatives among genotyped individuals
        int n = dosages.length;
        class Relative implements Comparable<Relative> {
            int idx; double sim;
            Relative(int idx, double sim) { this.idx = idx; this.sim = sim; }
            public int compareTo(Relative o) { return Double.compare(o.sim, this.sim); } // Descending
        }
        
        PriorityQueue<Relative> pq = new PriorityQueue<>();
        int originalTargetIdx = validIndices.get(targetIdx);

        for (int i = 0; i < n; i++) {
            if (Double.isNaN(dosages[i])) continue;
            int originalOtherIdx = validIndices.get(i);
            double sim = kinship[originalTargetIdx][originalOtherIdx];
            pq.add(new Relative(i, sim));
        }

        if (pq.isEmpty()) return ploidy / 2.0; // Fallback

        double weightedSum = 0;
        double totalWeight = 0;
        int count = 0;
        while (!pq.isEmpty() && count < k) {
            Relative r = pq.poll();
            double weight = Math.max(0.001, r.sim + 1.0); // Shift kinship to positive weight
            weightedSum += dosages[r.idx] * weight;
            totalWeight += weight;
            count++;
        }
        return weightedSum / totalWeight;
    }

    private double[] recode(double[] X, double threshold) {
        double[] res = new double[X.length];
        for (int i = 0; i < X.length; i++) res[i] = X[i] >= threshold ? 1.0 : 0.0;
        return res;
    }

    private double[] recodeRef(double[] X, double threshold) {
        double[] res = new double[X.length];
        for (int i = 0; i < X.length; i++) res[i] = X[i] <= threshold ? 1.0 : 0.0;
        return res;
    }

    public List<GwasInteraction> runEpistasisScan(List<GwasHit> topHits, Map<String, List<MarkerData>> chromToMarkers, double[] Yf, DMatrixRMaj W, double[] weights, DMatrixRMaj U, DMatrixRMaj Ystar, DMatrixRMaj Wstar) {
        System.out.println("[GWAS] Starting targeted epistasis scan (Leads vs Genome)...");
        List<GwasInteraction> interactions = Collections.synchronizedList(new ArrayList<>());
        
        // Take top 5 leads
        List<GwasHit> leads = topHits.stream().limit(5).filter(h -> h.pValue < 1e-4).toList();
        if (leads.isEmpty()) return interactions;

        // Find MarkerData for leads
        Map<String, double[]> leadDosages = new HashMap<>();
        for (GwasHit h : leads) {
            for (List<MarkerData> list : chromToMarkers.values()) {
                for (MarkerData md : list) {
                    if (md.id.equals(h.markerId)) leadDosages.put(md.id, md.dosages);
                }
            }
        }

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        int n = Yf.length;
        int q = Wstar.numCols;

        for (String leadId : leadDosages.keySet()) {
            double[] X1 = leadDosages.get(leadId);
            DMatrixRMaj X1star = new DMatrixRMaj(n, 1);
            CommonOps_DDRM.multTransA(U, new DMatrixRMaj(X1), X1star);

            for (List<MarkerData> markers : chromToMarkers.values()) {
                executor.submit(() -> {
                    // Pre-allocate matrices per thread
                    DMatrixRMaj M = new DMatrixRMaj(n, q + 3); // W + X1 + X2 + Xint
                    DMatrixRMaj LHS = new DMatrixRMaj(q + 3, q + 3);
                    DMatrixRMaj RHS = new DMatrixRMaj(q + 3, 1);
                    DMatrixRMaj beta = new DMatrixRMaj(q + 3, 1);
                    DMatrixRMaj invLHS = new DMatrixRMaj(q + 3, q + 3);
                    DMatrixRMaj X2mat = new DMatrixRMaj(n, 1);
                    DMatrixRMaj X2star = new DMatrixRMaj(n, 1);
                    DMatrixRMaj Xintmat = new DMatrixRMaj(n, 1);
                    DMatrixRMaj Xintstar = new DMatrixRMaj(n, 1);

                    for (MarkerData md : markers) {
                        if (md.id.equals(leadId)) continue;
                        
                        double[] X2 = md.dosages;
                        double[] Xint = new double[n];
                        for (int k = 0; k < n; k++) Xint[k] = X1[k] * X2[k];

                        CommonOps_DDRM.multTransA(U, new DMatrixRMaj(X2), X2star);
                        CommonOps_DDRM.multTransA(U, new DMatrixRMaj(Xint), Xintstar);

                        CommonOps_DDRM.insert(Wstar, M, 0, 0);
                        CommonOps_DDRM.insert(X1star, M, 0, q);
                        CommonOps_DDRM.insert(X2star, M, 0, q + 1);
                        CommonOps_DDRM.insert(Xintstar, M, 0, q + 2);

                        LHS.zero(); RHS.zero();
                        for (int k = 0; k < n; k++) {
                            double w = weights[k];
                            for (int r = 0; r < q + 3; r++) {
                                RHS.add(r, 0, M.get(k, r) * w * Ystar.get(k, 0));
                                for (int c = 0; c < q + 3; c++) LHS.add(r, c, M.get(k, r) * w * M.get(k, c));
                            }
                        }

                        if (CommonOps_DDRM.solve(LHS, RHS, beta)) {
                            CommonOps_DDRM.invert(LHS, invLHS);
                            double bInt = beta.get(q + 2, 0);
                            double rssFull = 0;
                            for (int k = 0; k < n; k++) {
                                double pred = 0;
                                for (int r = 0; r < q + 3; r++) pred += M.get(k, r) * beta.get(r, 0);
                                rssFull += Math.pow(Ystar.get(k, 0) - pred, 2) * weights[k];
                            }
                            double sigma2 = rssFull / (n - q - 3);
                            double se = Math.sqrt(Math.abs(invLHS.get(q + 2, q + 2) * sigma2));
                            double t = bInt / (se + 1e-15);
                            double p = 2.0 * (1.0 - GwasMathUtils.tCDF(Math.abs(t), n - q - 3));

                            if (p < 1e-4) {
                                GwasInteraction gi = new GwasInteraction();
                                gi.marker1 = leadId; gi.marker2 = md.id;
                                gi.pValue = p; gi.effect = bInt;
                                interactions.add(gi);
                            }
                        }
                    }
                });
            }
        }
        executor.shutdown();
        try { executor.awaitTermination(1, TimeUnit.HOURS); } catch (InterruptedException e) { e.printStackTrace(); }
        
        interactions.sort((a, b) -> Double.compare(a.pValue, b.pValue));
        return interactions;
    }

    private void applyFDRCorrection(List<GwasHit> hits) {
        int m = hits.size();
        for (int i = 0; i < m; i++) {
            double qVal = hits.get(i).pValue * m / (i + 1.0);
            hits.get(i).qValue = Math.min(1.0, qVal);
        }
        // Ensure monotonicity
        for (int i = m - 2; i >= 0; i--) {
            hits.get(i).qValue = Math.min(hits.get(i).qValue, hits.get(i + 1).qValue);
        }
    }

    private void selectBestModels(List<GwasHit> hits) {
        Map<String, GwasHit> bestByMarker = new HashMap<>();
        for (GwasHit hit : hits) {
            GwasHit currentBest = bestByMarker.get(hit.markerId);
            if (currentBest == null || hit.aic < currentBest.aic) {
                bestByMarker.put(hit.markerId, hit);
            }
        }
        for (GwasHit hit : hits) {
            if (hit == bestByMarker.get(hit.markerId)) {
                hit.isBestModel = true;
            }
        }
    }

    private void pruneByLD(Map<String, List<MarkerData>> chromToMarkers) {
        for (String chrom : chromToMarkers.keySet()) {
            List<MarkerData> markers = chromToMarkers.get(chrom);
            if (markers.size() < 2) continue;
            
            List<MarkerData> pruned = new ArrayList<>();
            pruned.add(markers.get(0));
            
            for (int i = 1; i < markers.size(); i++) {
                MarkerData current = markers.get(i);
                MarkerData last = pruned.get(pruned.size() - 1);
                
                double r2 = calculateR2(current.dosages, last.dosages);
                if (r2 < ldPruneThreshold) {
                    pruned.add(current);
                }
            }
            chromToMarkers.put(chrom, pruned);
        }
    }

    private double calculateR2(double[] x, double[] y) {
        double sumX = 0, sumY = 0, sumXX = 0, sumYY = 0, sumXY = 0;
        int n = x.length;
        for (int i = 0; i < n; i++) {
            sumX += x[i]; sumY += y[i];
            sumXX += x[i] * x[i]; sumYY += y[i] * y[i];
            sumXY += x[i] * y[i];
        }
        double cov = (n * sumXY - sumX * sumY);
        double varX = (n * sumXX - sumX * sumX);
        double varY = (n * sumYY - sumY * sumY);
        if (varX == 0 || varY == 0) return 0;
        return (cov * cov) / (varX * varY);
    }
}
