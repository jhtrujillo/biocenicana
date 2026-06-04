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
        ADDITIVE, SIMPLEX_DOMINANT, DUPLEX_DOMINANT, TRIPLEX_DOMINANT,
        SIMPLEX_DOMINANT_REF, DUPLEX_DOMINANT_REF, TRIPLEX_DOMINANT_REF, GENERAL,
        DIPLO_ADDITIVE, DIPLO_GENERAL
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
        public String model;
        public boolean isBestModel;
        public List<Double>[] phenotypesByDosage;
        public List<String>[] samplesByDosage;
    }

    public static class GwasResult {
        public List<GwasHit> hits;
        public List<GwasInteraction> interactions = new ArrayList<>();
        public int totalMarkersScanned;
        public List<double[]> ldDecayData; // [distance, r2]
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
        boolean passesGenoFreq = true;
    }

    private int ploidy;
    private String[] sampleNames;
    private double[][] covariates;
    private double[][] kinship;
    private boolean useLoco = false;
    private boolean runEpistasis = false;
    private int windowSize = 1;
    private int numThreads = Runtime.getRuntime().availableProcessors();
    private Set<GeneticModel> models = new HashSet<>(List.of(GeneticModel.ADDITIVE));
    private double ldPruneThreshold = 1.0;
    private double maxMissing = 0.1;
    private double mafThreshold = 0.01;
    private double maxGenoFreq = 1.0;
    private String imputationMode = "mean";
    private boolean gwasPolyCompatibility = false;
    // FIX: reset per run() call to avoid accumulation across multiple traits
    private int totalMarkersUnfiltered = 0;

    public GwasEngine(int ploidy, String[] sampleNames) {
        this.ploidy = ploidy;
        this.sampleNames = sampleNames;
    }

    public void setKinship(double[][] kinship) { this.kinship = kinship; }
    public void setLoco(boolean useLoco) { this.useLoco = useLoco; }
    public void setRunEpistasis(boolean runEpistasis) { this.runEpistasis = runEpistasis; }
    public void setWindowSize(int windowSize) { this.windowSize = windowSize; }
    public void setModels(Set<GeneticModel> models) { this.models = models; }
    public void setLdPruneThreshold(double ldPruneThreshold) { this.ldPruneThreshold = ldPruneThreshold; }
    public void setImputationMode(String mode) { this.imputationMode = mode; }
    public void setMaxMissing(double maxMissing) { this.maxMissing = maxMissing; }
    public void setMafThreshold(double mafThreshold) { this.mafThreshold = mafThreshold; }
    public void setMaxGenoFreq(double maxGenoFreq) { this.maxGenoFreq = maxGenoFreq; }
    public void setCovariates(double[][] covariates) { this.covariates = covariates; }
    public void setNumThreads(int numThreads) { this.numThreads = numThreads; }

    public void setGwasPolyCompatibility(boolean compatibility) {
        this.gwasPolyCompatibility = compatibility;
        if (compatibility) this.imputationMode = "mean";
    }

    public double[][] combineAllCovariates(double[][] pcaCovar, double[][] qCovar, PhenotypeData pheno, List<String> fixedCols) {
        int n = sampleNames.length;
        int pcaCount = (pcaCovar != null) ? pcaCovar[0].length : 0;
        int qCount = (qCovar != null) ? qCovar[0].length : 0;

        List<double[]> dummyEncoded = new ArrayList<>();
        if (fixedCols != null) {
            for (String colName : fixedCols) {
                Map<String, List<Integer>> levelIndices = new HashMap<>();
                for (int i = 0; i < n; i++) {
                    Map<String, String> fx = pheno.getFixedEffects(sampleNames[i]);
                    String val = (fx != null) ? fx.get(colName) : null;
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
        // FIX: reset counter so multiple run() calls on the same instance don't accumulate
        totalMarkersUnfiltered = 0;

        System.out.println("[GWAS] Loading variants and filtering...");

        // 1. Filter individuals with missing phenotypes
        List<Integer> validIndices = new ArrayList<>();
        for (int i = 0; i < yValues.length; i++) {
            if (!Double.isNaN(yValues[i])) validIndices.add(i);
        }
        int nFiltered = validIndices.size();
        double[] Yf = new double[nFiltered];
        for (int i = 0; i < nFiltered; i++) Yf[i] = yValues[validIndices.get(i)];

        // 2. Load and filter markers
        Map<String, List<MarkerData>> chromToMarkers = new LinkedHashMap<>();
        int totalLoaded = 0;
        List<double[]> allDosagesForKinship = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(vcfPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                totalMarkersUnfiltered++;
                String[] cols = line.split("\t");
                double[] dosages = extractDosages(cols, validIndices);
                if (dosages == null) continue;

                allDosagesForKinship.add(dosages);

                MarkerData md = new MarkerData();
                md.id = cols[2].equals(".") ? cols[0] + ":" + cols[1] : cols[2];
                md.chrom = cols[0];
                md.pos = Long.parseLong(cols[1]);
                md.ref = cols[3];
                md.alt = cols[4];
                md.dosages = dosages;
                md.passesGenoFreq = checkGenoFreq(dosages);

                chromToMarkers.computeIfAbsent(md.chrom, k -> new ArrayList<>()).add(md);
                totalLoaded++;
            }
        }

        // LD Pruning
        if (ldPruneThreshold < 0.999) {
            System.out.println("[GWAS] Pruning markers by LD (r2 < " + ldPruneThreshold + ")...");
            int before = totalLoaded;
            pruneByLD(chromToMarkers);
            totalLoaded = chromToMarkers.values().stream().mapToInt(List::size).sum();
            System.out.println("[GWAS] Pruning complete: " + before + " -> " + totalLoaded + " markers.");
        }

        // Smart grouping for LOCO
        Map<String, List<MarkerData>> locoGroups = new LinkedHashMap<>();
        if (useLoco) {
            List<MarkerData> smallContigs = new ArrayList<>();
            for (Map.Entry<String, List<MarkerData>> entry : chromToMarkers.entrySet()) {
                if (entry.getValue().size() < 100) smallContigs.addAll(entry.getValue());
                else locoGroups.put(entry.getKey(), entry.getValue());
            }
            if (!smallContigs.isEmpty()) locoGroups.put("__SMALL_CONTIGS__", smallContigs);
        } else {
            locoGroups = chromToMarkers;
        }

        System.out.println("[GWAS] Total markers: " + totalLoaded + " in " + chromToMarkers.size() + " sequences.");
        if (useLoco) System.out.println("[GWAS] LOCO Groups formed: " + locoGroups.size());

        // 3. Prepare covariates (W)
        int numCovs = (covariates != null ? covariates[0].length : 0);
        DMatrixRMaj W = new DMatrixRMaj(nFiltered, 1 + numCovs);
        for (int i = 0; i < nFiltered; i++) W.set(i, 0, 1.0);
        if (covariates != null) {
            for (int i = 0; i < nFiltered; i++) {
                int originalIdx = validIndices.get(i);
                for (int j = 0; j < numCovs; j++) W.set(i, j + 1, covariates[originalIdx][j]);
            }
        }

        // 4. Global kinship if not provided
        if (kinship == null) {
            System.out.println("[GWAS] Computing global " + (gwasPolyCompatibility ? "GWASpoly" : "VanRaden") + " kinship matrix...");
            if (gwasPolyCompatibility) {
                kinship = PopulationStructureAnalyzer.calculateGwasPolyKinship(allDosagesForKinship);
            } else {
                List<double[]> kinshipDosages = new ArrayList<>();
                for (List<MarkerData> markers : chromToMarkers.values())
                    for (MarkerData md : markers) kinshipDosages.add(md.dosages);
                kinship = PopulationStructureAnalyzer.calculateVanRadenKinship(kinshipDosages, ploidy);
            }
        }

        List<GwasHit> allHits = Collections.synchronizedList(new ArrayList<>());
        String[] filteredNames = new String[nFiltered];
        for (int i = 0; i < nFiltered; i++) filteredNames[i] = sampleNames[validIndices.get(i)];

        if (!useLoco) {
            System.out.println("[GWAS] Running Global Analysis...");
            runChromBlock(chromToMarkers.keySet(), chromToMarkers, kinship, nFiltered, filteredNames, Yf, W, allHits);
        } else {
            System.out.println("[GWAS] Running LOCO Analysis (Leave-One-Chromosome-Out)...");
            for (String groupName : locoGroups.keySet()) {
                System.out.print("[GWAS] Analyzing Group: " + groupName + " (LOCO) -> ");

                List<double[]> dosagesExcl = new ArrayList<>();
                for (Map.Entry<String, List<MarkerData>> entry : locoGroups.entrySet()) {
                    if (!entry.getKey().equals(groupName))
                        for (MarkerData md : entry.getValue()) dosagesExcl.add(md.dosages);
                }

                double[][] locoKinship;
                if (gwasPolyCompatibility) {
                    locoKinship = PopulationStructureAnalyzer.calculateGwasPolyKinship(dosagesExcl);
                } else {
                    locoKinship = PopulationStructureAnalyzer.calculateVanRadenKinship(dosagesExcl, ploidy);
                }

                if (locoKinship == null) {
                    System.out.println("Using global kinship (no other markers to exclude).");
                    locoKinship = kinship;
                } else {
                    System.out.println("Using specific LOCO kinship.");
                }

                Map<String, List<MarkerData>> markersInGroup = new HashMap<>();
                if (groupName.equals("__SMALL_CONTIGS__")) {
                    for (MarkerData md : locoGroups.get(groupName))
                        markersInGroup.computeIfAbsent(md.chrom, k -> new ArrayList<>()).add(md);
                } else {
                    markersInGroup.put(groupName, locoGroups.get(groupName));
                }

                runChromBlock(markersInGroup.keySet(), markersInGroup, locoKinship, nFiltered, filteredNames, Yf, W, allHits);
            }
        }

        selectBestModels(allHits);
        applyFDRCorrection(allHits);
        allHits.sort((a, b) -> Double.compare(a.pValue, b.pValue));

        GwasResult result = new GwasResult();
        result.hits = allHits;
        result.totalMarkersScanned = totalMarkersUnfiltered;

        if (runEpistasis) {
            System.out.println("[GWAS] Running epistasis phase...");
            List<double[]> allDosages = new ArrayList<>();
            for (List<MarkerData> list : chromToMarkers.values())
                for (MarkerData md : list) allDosages.add(md.dosages);

            double[][] globalK = gwasPolyCompatibility
                    ? PopulationStructureAnalyzer.calculateGwasPolyKinship(allDosages)
                    : PopulationStructureAnalyzer.calculateVanRadenKinship(allDosages, ploidy);

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

            // FIX: estimate REML delta instead of using hardcoded 1.0
            double delta = estimateDeltaREML(Ystar, Wstar, lambdas);
            System.out.printf("[GWAS] Epistasis REML Delta: %.4f%n", delta);

            double[] weights = new double[nFiltered];
            for (int i = 0; i < nFiltered; i++) weights[i] = 1.0 / (Math.max(1e-6, lambdas[i]) + delta);

            result.interactions = runEpistasisScan(allHits, chromToMarkers, Yf, W, weights, U, Ystar, Wstar);
        }

        result.ldDecayData = computeSampledLD(chromToMarkers, 2000);
        System.out.println("[GWAS] Calculated " + result.ldDecayData.size() + " LD data points for decay plot.");
        return result;
    }

    private List<double[]> computeSampledLD(Map<String, List<MarkerData>> chromToMarkers, int samples) {
        List<double[]> data = new ArrayList<>();
        Random rand = new Random(42);
        List<MarkerData> all = new ArrayList<>();
        for (List<MarkerData> list : chromToMarkers.values()) all.addAll(list);
        if (all.size() < 2) return data;

        // FIX: cap total attempts to avoid infinite loop when few same-chrom pairs exist
        int maxAttempts = samples * 50;
        int attempts = 0;
        while (data.size() < samples && attempts < maxAttempts) {
            attempts++;
            int i = rand.nextInt(all.size());
            int j = rand.nextInt(all.size());
            if (i == j) continue;
            MarkerData m1 = all.get(i);
            MarkerData m2 = all.get(j);
            if (!m1.chrom.equals(m2.chrom)) continue;
            long dist = Math.abs(m1.pos - m2.pos);
            if (dist > 30_000_000) continue;

            double r2 = GwasMathUtils.calculateR2(m1.dosages, m2.dosages);
            data.add(new double[]{dist, r2});
        }
        data.sort((a, b) -> Double.compare(a[0], b[0]));
        return data;
    }

    public int getTotalMarkersUnfiltered() {
        return totalMarkersUnfiltered;
    }

    private void runChromBlock(Set<String> chromsToScan, Map<String, List<MarkerData>> chromToMarkers,
                               double[][] kMatrix, int nFiltered, String[] filteredNames,
                               double[] Yf, DMatrixRMaj W, List<GwasHit> allHits) {
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

        DMatrixRMaj Ymat = new DMatrixRMaj(nFiltered, 1);
        for (int i = 0; i < nFiltered; i++) Ymat.set(i, 0, Yf[i]);
        DMatrixRMaj Ystar = new DMatrixRMaj(nFiltered, 1);
        CommonOps_DDRM.multTransA(U, Ymat, Ystar);
        DMatrixRMaj Wstar = new DMatrixRMaj(nFiltered, W.numCols);
        CommonOps_DDRM.multTransA(U, W, Wstar);

        double delta = estimateDeltaREML(Ystar, Wstar, lambdas);
        System.out.printf("[GWAS] Estimated REML Delta (VarE / VarG): %.4f%n", delta);

        double[] weights = new double[nFiltered];
        for (int i = 0; i < nFiltered; i++) weights[i] = 1.0 / (Math.max(1e-6, lambdas[i]) + delta);
        final double rssReduced = calculateRSS(Ystar, Wstar, lambdas, delta);

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        for (String chrom : chromsToScan) {
            List<MarkerData> markers = chromToMarkers.get(chrom);
            if (markers == null) continue;

            for (int i = 0; i < markers.size(); i += windowSize) {
                final List<MarkerData> window = markers.subList(i, Math.min(i + windowSize, markers.size()));
                executor.submit(() -> {
                    if (window.size() == 1) {
                        MarkerData md = window.get(0);
                        if (!md.passesGenoFreq) return;
                        for (GeneticModel model : models) {
                            double[] X = recodeForModel(md.dosages, model);
                            if (X == null) continue;
                            int df = (model == GeneticModel.GENERAL) ? getUniqueDosageCount(md.dosages) - 1 : 1;
                            GwasHit hit = testModel(md.id, md.chrom, md.pos, md.ref, md.alt,
                                    X, df, filteredNames, Yf, Ystar, Wstar, U, weights, rssReduced, model.toString());
                            if (hit.pValue < 1.0) allHits.add(hit);
                        }
                    } else {
                        double[] haploDosage = computeHaploDosage(window);
                        MarkerData first = window.get(0);
                        MarkerData last = window.get(window.size() - 1);
                        String blockId = String.format("Block_%s:%d-%d", first.chrom, first.pos, last.pos);
                        GwasHit hit = testModel(blockId, first.chrom, (first.pos + last.pos) / 2,
                                "Block", "Block", haploDosage, 1, filteredNames, Yf, Ystar, Wstar,
                                U, weights, rssReduced, "Haplotype");
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

    private double estimateDeltaREML(DMatrixRMaj Ystar, DMatrixRMaj Wstar, double[] lambdas) {
        int n = Ystar.numRows;
        int q = Wstar.numCols;

        java.util.function.DoubleFunction<Double> negLogLikelihood = (delta) -> {
            double rss = calculateRSS(Ystar, Wstar, lambdas, delta);
            if (rss >= 1e10 || rss <= 0) return 1e10;
            double sumLogH = 0;
            DMatrixRMaj WtDW = new DMatrixRMaj(q, q);
            for (int i = 0; i < n; i++) {
                double w = 1.0 / (lambdas[i] + delta);
                sumLogH += Math.log(lambdas[i] + delta);
                for (int j = 0; j < q; j++)
                    for (int k = 0; k < q; k++)
                        WtDW.add(j, k, Wstar.get(i, j) * w * Wstar.get(i, k));
            }
            double detWtDW = CommonOps_DDRM.det(WtDW);
            if (detWtDW <= 0) return 1e10;
            return (n - q) * Math.log(rss) + sumLogH + Math.log(detWtDW);
        };

        return GwasMathUtils.brentOptimize(1e-4, 1000.0, 1e-4, negLogLikelihood);
    }

    // FIX: removed 6 unused pre-allocated matrix parameters (M, MtVinv, LHS, RHS, beta, invLHS)
    private GwasHit testModel(String id, String chrom, long pos, String ref, String alt,
                              double[] X, int df, String[] filteredNames, double[] Yf,
                              DMatrixRMaj Ystar, DMatrixRMaj Wstar, DMatrixRMaj U,
                              double[] weights, double rssReduced, String modelName) {
        int n = Yf.length;
        int q = Wstar.numCols;
        int markerCols = df;

        DMatrixRMaj M_multi = new DMatrixRMaj(n, q + markerCols);
        CommonOps_DDRM.insert(Wstar, M_multi, 0, 0);

        DMatrixRMaj Xmat = new DMatrixRMaj(n, 1);
        DMatrixRMaj Xstar = new DMatrixRMaj(n, 1);
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
        double fStat = ((rssReduced - rssFull) / df) / (rssFull / (n - q - df));
        double p = 1.0 - GwasMathUtils.fCDF(fStat, df, n - q - df);
        double aic = n * Math.log(rssFull / n) + 2 * (q + df + 1);

        GwasHit hit = new GwasHit();
        hit.markerId = id; hit.chromosome = chrom; hit.position = pos;
        hit.refAllele = ref; hit.altAllele = alt;
        hit.effect = beta_multi.get(q, 0);
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
                int dosage = (int) Math.round(X[i * markerCols]);
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
            case TRIPLEX_DOMINANT: return recode(dosages, 3.0);
            case SIMPLEX_DOMINANT_REF: return recodeRef(dosages, ploidy - 1);
            case DUPLEX_DOMINANT_REF: return recodeRef(dosages, ploidy - 2);
            case TRIPLEX_DOMINANT_REF: return recodeRef(dosages, ploidy - 3);
            case GENERAL: return recodeGeneral(dosages);
            case DIPLO_ADDITIVE: {
                double[] da = new double[n];
                for (int i = 0; i < n; i++) {
                    double d = Math.round(dosages[i]);
                    da[i] = (d > 0 && d < ploidy) ? 1.0 : (d == ploidy ? 2.0 : 0.0);
                }
                return da;
            }
            case DIPLO_GENERAL: {
                double[] dg = new double[n];
                for (int i = 0; i < n; i++) {
                    double d = Math.round(dosages[i]);
                    dg[i] = (d > 0 && d < ploidy) ? 1.0 : (d == ploidy ? 2.0 : 0.0);
                }
                return recodeGeneral(dg);
            }
            default: return dosages;
        }
    }

    private double[] recodeGeneral(double[] dosages) {
        int n = dosages.length;
        Set<Integer> unique = new TreeSet<>();
        for (double d : dosages) unique.add((int) Math.round(d));
        if (unique.size() < 2) return null;
        List<Integer> levels = new ArrayList<>(unique);
        int df = levels.size() - 1;
        double[] Xmulti = new double[n * df];
        for (int i = 0; i < n; i++) {
            int d = (int) Math.round(dosages[i]);
            for (int j = 0; j < df; j++) Xmulti[i * df + j] = (d == levels.get(j)) ? 1.0 : 0.0;
        }
        return Xmulti;
    }

    private int getUniqueDosageCount(double[] dosages) {
        Set<Integer> unique = new HashSet<>();
        for (double d : dosages) unique.add((int) Math.round(d));
        return unique.size();
    }

    /**
     * Extracts and imputes dosage values for valid samples from a VCF data row.
     * Returns null if the marker fails missing-rate or MAF filters.
     */
    private double[] extractDosages(String[] cols, List<Integer> validIndices) {
        int nFiltered = validIndices.size();
        double[] dosages = new double[nFiltered];
        List<Integer> missingIndices = new ArrayList<>();
        double sum = 0;
        int count = 0;

        // FIX: skip multiallelic in gwasPolyCompatibility mode
        if (gwasPolyCompatibility && cols[4].contains(",")) return null;

        // Resolve DS field index once for this variant
        int dsIndex = -1;
        if (cols.length > 8) {
            String[] formatFields = cols[8].split(":");
            for (int f = 0; f < formatFields.length; f++) {
                if (formatFields[f].equalsIgnoreCase("DS")) { dsIndex = f; break; }
            }
        }

        for (int i = 0; i < nFiltered; i++) {
            int idx = validIndices.get(i);
            String gData = cols[9 + idx];
            if (gData.startsWith(".")) {
                missingIndices.add(i);
                dosages[i] = Double.NaN;
                continue;
            }

            String[] fields = gData.split(":");
            double dosage;
            if (dsIndex != -1 && fields.length > dsIndex) {
                try {
                    dosage = Double.parseDouble(fields[dsIndex]);
                } catch (NumberFormatException e) {
                    // FIX: only count '1' alleles — correct for biallelic; avoids wrong counts for '2', '3', etc.
                    dosage = countAltAlleles(fields[0]);
                }
            } else {
                dosage = countAltAlleles(fields[0]);
            }
            dosages[i] = dosage;
            sum += dosage;
            count++;
        }

        double missingRatio = (double) missingIndices.size() / nFiltered;
        if (count == 0 || missingRatio > maxMissing) return null;
        if (maxMissing == 0 && !missingIndices.isEmpty()) return null;

        double mean = sum / count;
        double freq = sum / (count * ploidy);
        double mafValue = Math.min(freq, 1.0 - freq);
        if (mafValue < mafThreshold) return null;

        // Impute missing values
        for (int i : missingIndices) {
            if ("knn".equalsIgnoreCase(imputationMode) && kinship != null) {
                dosages[i] = imputeWithKNN(i, dosages, validIndices, kinship, 5);
            } else {
                dosages[i] = mean;
            }
        }
        return dosages;
    }

    // FIX: count only '1' alleles — safe for biallelic VCF; avoids miscounting multiallelic alleles (0/2 → 0, not 1)
    private int countAltAlleles(String gt) {
        int count = 0;
        for (char c : gt.toCharArray()) if (c == '1') count++;
        return count;
    }

    private boolean checkGenoFreq(double[] dosages) {
        if (maxGenoFreq >= 1.0) return true;
        Map<Integer, Integer> counts = new HashMap<>();
        int total = 0;
        for (double d : dosages) {
            if (!Double.isNaN(d)) {
                int id = (int) Math.round(d);
                counts.put(id, counts.getOrDefault(id, 0) + 1);
                total++;
            }
        }
        if (total == 0) return false;
        for (int c : counts.values()) if ((double) c / total > maxGenoFreq) return false;
        return true;
    }

    private double imputeWithKNN(int targetIdx, double[] dosages, List<Integer> validIndices, double[][] kinship, int k) {
        int n = dosages.length;
        class Relative implements Comparable<Relative> {
            int idx; double sim;
            Relative(int idx, double sim) { this.idx = idx; this.sim = sim; }
            public int compareTo(Relative o) { return Double.compare(o.sim, this.sim); }
        }
        PriorityQueue<Relative> pq = new PriorityQueue<>();
        int originalTargetIdx = validIndices.get(targetIdx);
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(dosages[i])) continue;
            int originalOtherIdx = validIndices.get(i);
            pq.add(new Relative(i, kinship[originalTargetIdx][originalOtherIdx]));
        }
        if (pq.isEmpty()) return ploidy / 2.0;
        double weightedSum = 0, totalWeight = 0;
        int count = 0;
        while (!pq.isEmpty() && count < k) {
            Relative r = pq.poll();
            double weight = Math.max(0.001, r.sim + 1.0);
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

    public List<GwasInteraction> runEpistasisScan(List<GwasHit> topHits, Map<String, List<MarkerData>> chromToMarkers,
                                                   double[] Yf, DMatrixRMaj W, double[] weights, DMatrixRMaj U,
                                                   DMatrixRMaj Ystar, DMatrixRMaj Wstar) {
        System.out.println("[GWAS] Starting targeted epistasis scan (Leads vs Genome)...");
        List<GwasInteraction> interactions = Collections.synchronizedList(new ArrayList<>());

        List<GwasHit> leads = topHits.stream().limit(5).filter(h -> h.pValue < 1e-4).toList();
        if (leads.isEmpty()) return interactions;

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
                    DMatrixRMaj M = new DMatrixRMaj(n, q + 3);
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

    private void selectBestModels(List<GwasHit> hits) {
        for (GwasHit h : hits) h.isBestModel = false;
        Map<String, GwasHit> bestByMarker = new HashMap<>();
        for (GwasHit hit : hits) {
            GwasHit currentBest = bestByMarker.get(hit.markerId);
            if (currentBest == null || hit.pValue < currentBest.pValue) bestByMarker.put(hit.markerId, hit);
        }
        for (GwasHit hit : hits) {
            if (hit == bestByMarker.get(hit.markerId)) hit.isBestModel = true;
        }
    }

    private void applyFDRCorrection(List<GwasHit> hits) {
        List<GwasHit> bestHits = new ArrayList<>();
        for (GwasHit h : hits) if (h.isBestModel) bestHits.add(h);
        bestHits.sort((a, b) -> Double.compare(a.pValue, b.pValue));
        int m = bestHits.size();
        for (int i = 0; i < m; i++) {
            bestHits.get(i).qValue = Math.min(1.0, bestHits.get(i).pValue * m / (i + 1.0));
        }
        for (int i = m - 2; i >= 0; i--) {
            bestHits.get(i).qValue = Math.min(bestHits.get(i).qValue, bestHits.get(i + 1).qValue);
        }
        for (GwasHit h : hits) if (!h.isBestModel) h.qValue = 1.0;
    }

    /**
     * LD pruning using a sliding window: a candidate is retained only if its r² with
     * ALL previously retained markers in the window is below the threshold.
     */
    private void pruneByLD(Map<String, List<MarkerData>> chromToMarkers) {
        for (String chrom : chromToMarkers.keySet()) {
            List<MarkerData> markers = chromToMarkers.get(chrom);
            if (markers.size() < 2) continue;

            List<MarkerData> pruned = new ArrayList<>();
            pruned.add(markers.get(0));

            for (int i = 1; i < markers.size(); i++) {
                MarkerData current = markers.get(i);
                boolean keep = true;
                // FIX: compare against all retained markers, not just the last one
                for (MarkerData retained : pruned) {
                    if (GwasMathUtils.calculateR2(current.dosages, retained.dosages) >= ldPruneThreshold) {
                        keep = false;
                        break;
                    }
                }
                if (keep) pruned.add(current);
            }
            chromToMarkers.put(chrom, pruned);
        }
    }
}
