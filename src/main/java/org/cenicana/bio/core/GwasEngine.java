package org.cenicana.bio.core;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
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

    public static class GwasHit {
        public String markerId;
        public String chromosome;
        public long position;
        public String refAllele;
        public String altAllele;
        public double pValue;
        public double effect;
        public double r2;
        public String model; // Additive, SimplexDominant, etc.
        public List<Double>[] phenotypesByDosage; // To store distributions for boxplots
        public List<String>[] samplesByDosage; // To track elite candidates
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
    private int numThreads = Runtime.getRuntime().availableProcessors();

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

    public void setCovariates(double[][] covariates) { this.covariates = covariates; }
    public void setNumThreads(int numThreads) { this.numThreads = numThreads; }

    public double[][] combineCovariates(double[][] pcaCovar, PhenotypeData pheno, List<String> fixedCols) {
        if (fixedCols == null || fixedCols.isEmpty()) return pcaCovar;
        
        List<double[]> dummyEncoded = new ArrayList<>();
        
        for (String colName : fixedCols) {
            Map<String, List<Integer>> levelIndices = new HashMap<>();
            
            for (int i = 0; i < sampleNames.length; i++) {
                String val = pheno.getFixedEffects(sampleNames[i]).get(colName);
                if (val != null) {
                    levelIndices.computeIfAbsent(val, k -> new ArrayList<>()).add(i);
                }
            }
            
            if (levelIndices.isEmpty()) continue;
            
            // Create n-1 columns for this factor
            List<String> uniqueLevels = new ArrayList<>(levelIndices.keySet());
            for (int i = 0; i < uniqueLevels.size() - 1; i++) {
                double[] dummy = new double[sampleNames.length];
                for (int idx : levelIndices.get(uniqueLevels.get(i))) dummy[idx] = 1.0;
                dummyEncoded.add(dummy);
            }
        }
        
        int n = sampleNames.length;
        int pcaCount = (pcaCovar != null) ? pcaCovar[0].length : 0;
        int totalCovs = pcaCount + dummyEncoded.size();
        
        double[][] combined = new double[n][totalCovs];
        for (int i = 0; i < n; i++) {
            if (pcaCovar != null) System.arraycopy(pcaCovar[i], 0, combined[i], 0, pcaCount);
            for (int j = 0; j < dummyEncoded.size(); j++) {
                combined[i][pcaCount + j] = dummyEncoded.get(j)[i];
            }
        }
        return combined;
    }

    public List<GwasHit> run(String vcfPath, double[] yValues, String traitName) throws Exception {
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
                double[] dosages = extractDosages(cols, validIndices);
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
        return allHits;
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
            
            int batchSize = 1000;
            for (int i = 0; i < markers.size(); i += batchSize) {
                final List<MarkerData> batch = markers.subList(i, Math.min(i + batchSize, markers.size()));
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

                    for (MarkerData md : batch) {
                        GwasHit hit = testModel(md.id, md.chrom, md.pos, md.ref, md.alt, md.dosages, filteredNames, Yf, Ystar, Wstar, U, weights, rssReduced, "Additive", M, MtVinv, LHS, RHS, beta, invLHS, Xmat, Xstar);
                        allHits.add(hit);
                    }
                });
            }
        }
        executor.shutdown();
        try { executor.awaitTermination(1, TimeUnit.HOURS); } catch (InterruptedException e) { e.printStackTrace(); }
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

    private GwasHit testModel(String id, String chrom, long pos, String ref, String alt, double[] X, String[] filteredNames, double[] Yf, DMatrixRMaj Ystar, DMatrixRMaj Wstar, DMatrixRMaj U, double[] weights, double rssReduced, String modelName,
                              DMatrixRMaj M, DMatrixRMaj MtVinv, DMatrixRMaj LHS, DMatrixRMaj RHS, DMatrixRMaj beta, DMatrixRMaj invLHS, DMatrixRMaj Xmat, DMatrixRMaj Xstar) {
        int n = X.length;
        int q = Wstar.numCols;
        for(int i=0; i<n; i++) Xmat.set(i, 0, X[i]);
        CommonOps_DDRM.multTransA(U, Xmat, Xstar);
        CommonOps_DDRM.insert(Wstar, M, 0, 0);
        CommonOps_DDRM.insert(Xstar, M, 0, q);

        LHS.zero(); RHS.zero();
        for (int i = 0; i < n; i++) {
            double w = weights[i];
            for (int j = 0; j < q + 1; j++) {
                double val = M.get(i, j) * w;
                RHS.add(j, 0, val * Ystar.get(i, 0));
                for (int k = 0; k < q + 1; k++) LHS.add(j, k, val * M.get(i, k));
            }
        }

        if (!CommonOps_DDRM.solve(LHS, RHS, beta)) {
            GwasHit h = new GwasHit(); h.pValue = 1.0; return h;
        }

        CommonOps_DDRM.invert(LHS, invLHS);
        double bSNP = beta.get(q, 0);
        double rssFull = 0;
        for (int i = 0; i < n; i++) {
            double pred = 0;
            for (int j = 0; j < q + 1; j++) pred += M.get(i, j) * beta.get(j, 0);
            double err = Ystar.get(i, 0) - pred;
            rssFull += (err * err) * weights[i];
        }

        double partialR2 = (rssReduced - rssFull) / (rssReduced + 1e-10);
        double sigma2g = rssFull / (n - q - 1);
        double se = Math.sqrt(Math.abs(invLHS.get(q, q) * sigma2g));
        double tStat = bSNP / (se + 1e-15);
        
        GwasHit hit = new GwasHit();
        hit.markerId = id; hit.chromosome = chrom; hit.position = pos;
        hit.refAllele = ref; hit.altAllele = alt;
        hit.effect = bSNP; hit.model = modelName; hit.r2 = Math.max(0, partialR2);
        double p = 2.0 * (1.0 - GwasMathUtils.tCDF(Math.abs(tStat), n - q - 1));
        hit.pValue = Math.max(p, 1e-300);
        
        if (hit.pValue < 0.001) {
            hit.phenotypesByDosage = new List[ploidy + 1];
            hit.samplesByDosage = new List[ploidy + 1];
            for (int i = 0; i <= ploidy; i++) {
                hit.phenotypesByDosage[i] = new ArrayList<>();
                hit.samplesByDosage[i] = new ArrayList<>();
            }
            for (int i = 0; i < n; i++) {
                int dosage = (int) Math.round(X[i]);
                if (dosage >= 0 && dosage <= ploidy) {
                    hit.phenotypesByDosage[dosage].add(Yf[i]);
                    hit.samplesByDosage[dosage].add(filteredNames[i]);
                }
            }
        }
        return hit;
    }

    private double[] extractDosages(String[] cols, List<Integer> validIndices) {
        int n = sampleNames.length;
        double[] dosages = new double[validIndices.size()];
        int missing = 0;
        double sum = 0;

        for (int i = 0; i < validIndices.size(); i++) {
            int idx = validIndices.get(i);
            String gData = cols[9 + idx];
            if (gData.startsWith(".")) {
                missing++;
                dosages[i] = Double.NaN;
            } else {
                String gt = gData.split(":")[0];
                int count = 0;
                for (char c : gt.toCharArray()) if (c > '0' && c <= '9') count++;
                dosages[i] = count;
                sum += count;
            }
        }

        double maf = sum / (validIndices.size() * ploidy);
        if (maf > 0.5) maf = 1.0 - maf;
        if (maf < 0.01 || (double)missing/n > 0.2) return null;

        // Simple imputation with mean
        double mean = sum / (validIndices.size() - missing);
        for (int i = 0; i < dosages.length; i++) {
            if (Double.isNaN(dosages[i])) dosages[i] = mean;
        }
        return dosages;
    }

    private double[] recode(double[] X, double threshold) {
        double[] res = new double[X.length];
        for (int i = 0; i < X.length; i++) res[i] = X[i] >= threshold ? 1.0 : 0.0;
        return res;
    }
}
