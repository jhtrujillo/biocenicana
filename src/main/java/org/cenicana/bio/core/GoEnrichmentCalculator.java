package org.cenicana.bio.core;

import java.util.*;

/**
 * Calculates Gene Ontology (GO) enrichment using a hypergeometric test.
 *
 * <p><b>Numerical Stability:</b> All combinatorial calculations are performed
 * in log-space using the Lanczos approximation of the log-gamma function,
 * avoiding floating-point overflow that occurs with large genomes (N > 10,000).
 *
 * <p><b>Multiple Testing Correction:</b> Applies Benjamini-Hochberg (BH)
 * False Discovery Rate correction, which is the standard in bioinformatics GO
 * enrichment analysis.
 */
public class GoEnrichmentCalculator {

    // Lanczos approximation coefficients (g=7, n=9)
    private static final double[] LANCZOS_COEFFICIENTS = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    public static class EnrichmentResult {
        public String goId;
        public String description;
        public int studyCount;
        public int backgroundCount;
        /** Raw hypergeometric p-value (P(X >= k)). */
        public double pValue;
        /** Benjamini-Hochberg FDR corrected p-value. */
        public double fdrCorrectedPValue;

        public EnrichmentResult(String id, int sc, int bc, double p) {
            this.goId = id;
            this.studyCount = sc;
            this.backgroundCount = bc;
            this.pValue = p;
            this.fdrCorrectedPValue = p; // Will be overwritten by BH correction
        }
    }

    /**
     * Runs GO enrichment analysis.
     *
     * @param studyGenes     Map of GeneID → List of GO terms for genes in the study set (e.g., a syntenic block)
     * @param backgroundGenes Map of GeneID → List of GO terms for all genes in the genome (background)
     * @return Sorted list of enriched GO terms with BH-corrected p-values
     */
    public List<EnrichmentResult> calculate(Map<String, List<String>> studyGenes,
                                            Map<String, List<String>> backgroundGenes) {
        int n = studyGenes.size();  // Study set size
        int N = backgroundGenes.size(); // Background size

        if (n == 0 || N == 0) return Collections.emptyList();

        Map<String, Integer> studyCounts = countGo(studyGenes);
        Map<String, Integer> backgroundCounts = countGo(backgroundGenes);

        // Step 1: Calculate raw hypergeometric p-values in parallel
        List<EnrichmentResult> results = Collections.synchronizedList(new ArrayList<>());
        studyCounts.keySet().parallelStream().forEach(goId -> {
            int k = studyCounts.get(goId);           // Successes in study
            int K = backgroundCounts.getOrDefault(goId, k); // Successes in background
            double pValue = hypergeometricPValue(N, K, n, k);
            if (pValue < 0.05) {
                results.add(new EnrichmentResult(goId, k, K, pValue));
            }
        });

        // Step 2: Sort by raw p-value (required before BH correction)
        results.sort(Comparator.comparingDouble(r -> r.pValue));

        // Step 3: Apply Benjamini-Hochberg FDR correction
        applyBenjaminiHochberg(results);

        return results;
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Counts the number of genes annotated with each GO term.
     * Each gene is counted only once per GO term (de-duplicated).
     */
    private Map<String, Integer> countGo(Map<String, List<String>> geneMap) {
        Map<String, Integer> counts = new HashMap<>();
        for (List<String> terms : geneMap.values()) {
            Set<String> uniqueTerms = new HashSet<>(terms);
            for (String term : uniqueTerms) {
                counts.merge(term, 1, Integer::sum);
            }
        }
        return counts;
    }

    /**
     * Calculates the one-tailed hypergeometric p-value P(X >= k) in log-space
     * to prevent numerical overflow for large N (e.g., sugarcane genomes with
     * ~35,000 genes across both haplomes).
     *
     * <p>Uses the log-sum-exp trick: sum terms in log-space and convert back at
     * the very end, maintaining full double precision throughout.
     *
     * @param N total population size (background genes)
     * @param K successes in population (background genes with this GO term)
     * @param n sample size (study genes in the block)
     * @param k observed successes (study genes with this GO term)
     * @return P(X >= k) under the hypergeometric distribution
     */
    private double hypergeometricPValue(int N, int K, int n, int k) {
        if (k <= 0) return 1.0;
        if (K > N || n > N) return 1.0;

        int limit = Math.min(n, K);

        // Compute each term log(P(X = i)) = logC(K,i) + logC(N-K, n-i) - logC(N,n)
        double logDenominator = logCombination(N, n);
        double[] logTerms = new double[limit - k + 1];
        for (int i = k; i <= limit; i++) {
            logTerms[i - k] = logCombination(K, i) + logCombination(N - K, n - i) - logDenominator;
        }

        return Math.min(1.0, Math.exp(logSumExp(logTerms)));
    }

    /**
     * Computes log(C(n, k)) = log(n!) - log(k!) - log((n-k)!) using the
     * Lanczos log-gamma approximation, safe for arbitrarily large n.
     */
    private double logCombination(int n, int k) {
        if (k < 0 || k > n) return Double.NEGATIVE_INFINITY;
        if (k == 0 || k == n) return 0.0;
        if (k > n - k) k = n - k; // Use symmetry to reduce iterations
        return logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1);
    }

    /**
     * Log-Gamma function using the Lanczos approximation.
     * Accurate to ~15 significant figures for x > 0.5.
     */
    private double logGamma(double x) {
        if (x < 0.5) {
            // Reflection formula: Γ(x)Γ(1-x) = π/sin(πx)
            return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1.0 - x);
        }
        x -= 1.0;
        double a = LANCZOS_COEFFICIENTS[0];
        double t = x + 7.5; // g = 7
        for (int i = 1; i < LANCZOS_COEFFICIENTS.length; i++) {
            a += LANCZOS_COEFFICIENTS[i] / (x + i);
        }
        return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
    }

    /**
     * Numerically stable log-sum-exp: log(Σ exp(terms)).
     * Subtracts the maximum term before summing to avoid underflow.
     */
    private double logSumExp(double[] logTerms) {
        if (logTerms.length == 0) return Double.NEGATIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (double v : logTerms) if (v > max) max = v;
        if (Double.isInfinite(max)) return max;
        double sum = 0.0;
        for (double v : logTerms) sum += Math.exp(v - max);
        return max + Math.log(sum);
    }

    /**
     * Applies Benjamini-Hochberg FDR correction in-place.
     * Assumes the list is already sorted by ascending raw p-value.
     *
     * <p>BH adjusted p-value for rank i (1-indexed) of m total tests:
     * <pre>p_adj(i) = p(i) * m / i</pre>
     * Applied in reverse (from largest rank) to ensure monotonicity.
     */
    private void applyBenjaminiHochberg(List<EnrichmentResult> results) {
        int m = results.size();
        if (m == 0) return;

        // Traverse from last (largest p) to first (smallest p)
        double runningMin = 1.0;
        for (int i = m - 1; i >= 0; i--) {
            double adjusted = results.get(i).pValue * m / (i + 1);
            runningMin = Math.min(runningMin, adjusted);
            results.get(i).fdrCorrectedPValue = runningMin;
        }
    }
}
