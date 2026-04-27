package org.cenicana.bio.core;

import java.util.*;

/**
 * Calculates Gene Ontology (GO) enrichment using a hypergeometric-like approach (simplified Fisher's Exact).
 * Given a set of genes in a block (study set) and all genes in the genome (background set).
 */
public class GoEnrichmentCalculator {

    public static class EnrichmentResult {
        public String goId;
        public String description;
        public int studyCount;
        public int backgroundCount;
        public double pValue;

        public EnrichmentResult(String id, int sc, int bc, double p) {
            this.goId = id; this.studyCount = sc; this.backgroundCount = bc; this.pValue = p;
        }
    }

    /**
     * @param studyGenes Map of GeneID -> List of GO terms for genes in the block
     * @param backgroundGenes Map of GeneID -> List of GO terms for all genes
     * @return List of enriched GO terms with p < 0.05
     */
    public List<EnrichmentResult> calculate(Map<String, List<String>> studyGenes, Map<String, List<String>> backgroundGenes) {
        int n = studyGenes.size(); // Total genes in study (block)
        int N = backgroundGenes.size(); // Total genes in background (genome)
        
        if (n == 0 || N == 0) return Collections.emptyList();

        Map<String, Integer> studyCounts = countGo(studyGenes);
        Map<String, Integer> backgroundCounts = countGo(backgroundGenes);

        // Parallel processing of GO terms
        List<EnrichmentResult> results = java.util.Collections.synchronizedList(new ArrayList<>());
        studyCounts.keySet().parallelStream().forEach(goId -> {
            int k = studyCounts.get(goId);
            int K = backgroundCounts.getOrDefault(goId, k);
            double pValue = calculateHypergeometricPValue(N, K, n, k);
            if (pValue < 0.05) {
                results.add(new EnrichmentResult(goId, k, K, pValue));
            }
        });

        results.sort(Comparator.comparingDouble(r -> r.pValue));
        return results;
    }

    private Map<String, Integer> countGo(Map<String, List<String>> geneMap) {
        Map<String, Integer> counts = new HashMap<>();
        for (List<String> terms : geneMap.values()) {
            Set<String> uniqueTerms = new HashSet<>(terms); // Avoid double counting per gene
            for (String term : uniqueTerms) {
                counts.put(term, counts.getOrDefault(term, 0) + 1);
            }
        }
        return counts;
    }

    /**
     * Simplified Hypergeometric P-value calculation (P(X >= k)).
     * N: total population, K: success in population, n: sample size, k: success in sample.
     */
    private double calculateHypergeometricPValue(int N, int K, int n, int k) {
        double p = 0;
        for (int i = k; i <= Math.min(n, K); i++) {
            p += combination(K, i) * combination(N - K, n - i) / combination(N, n);
            if (p > 1) return 1.0; // Numerical stability fallback
        }
        return Math.min(1.0, p);
    }

    private double combination(int n, int k) {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n / 2) k = n - k;
        
        double res = 1;
        for (int i = 1; i <= k; i++) {
            res *= (double) (n - i + 1) / i;
        }
        return res;
    }
}
