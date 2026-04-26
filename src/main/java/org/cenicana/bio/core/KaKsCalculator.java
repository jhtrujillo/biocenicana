package org.cenicana.bio.core;

import java.util.HashMap;
import java.util.Map;

/**
 * Native implementation of Ka/Ks calculation using Nei-Gojobori (1986) method.
 * Includes a basic Needleman-Wunsch aligner for CDS sequences.
 */
public class KaKsCalculator {

    private static final Map<String, String> GENETIC_CODE = new HashMap<>();
    static {
        String[][] table = {
            {"TTT","F"},{"TTC","F"},{"TTA","L"},{"TTG","L"},
            {"TCT","S"},{"TCC","S"},{"TCA","S"},{"TCG","S"},
            {"TAT","Y"},{"TAC","Y"},{"TAA","*"},{"TAG","*"},
            {"TGT","C"},{"TGC","C"},{"TGA","*"},{"TGG","W"},
            {"CTT","L"},{"CTC","L"},{"CTA","L"},{"CTG","L"},
            {"CCT","P"},{"CCC","P"},{"CCA","P"},{"CCG","P"},
            {"CAT","H"},{"CAC","H"},{"CAA","Q"},{"CAG","Q"},
            {"CGT","R"},{"CGC","R"},{"CGA","R"},{"CGG","R"},
            {"ATT","I"},{"ATC","I"},{"ATA","I"},{"ATG","M"},
            {"ACT","T"},{"ACC","T"},{"ACA","T"},{"ACG","T"},
            {"AAT","N"},{"AAC","N"},{"AAA","K"},{"AAG","K"},
            {"AGT","S"},{"AGC","S"},{"AGA","R"},{"AGG","R"},
            {"GTT","V"},{"GTC","V"},{"GTA","V"},{"GTG","V"},
            {"GCT","A"},{"GCC","A"},{"GCA","A"},{"GCG","A"},
            {"GAT","D"},{"GAC","D"},{"GAA","E"},{"GAG","E"},
            {"GGT","G"},{"GGC","G"},{"GGA","G"},{"GGG","G"}
        };
        for (String[] pair : table) GENETIC_CODE.put(pair[0], pair[1]);
    }

    public static class Result {
        public double ka, ks, ratio;
        public Result(double ka, double ks) {
            this.ka = ka;
            this.ks = ks;
            this.ratio = (ks > 0) ? ka / ks : 0;
        }
    }

    /**
     * Calculates Ka/Ks for two CDS sequences.
     * Sequences are assumed to be aligned or of same length.
     */
    public Result calculate(String seq1, String seq2) {
        seq1 = seq1.toUpperCase();
        seq2 = seq2.toUpperCase();
        
        // Simple validation: must be multiples of 3
        int len = Math.min(seq1.length(), seq2.length());
        len = (len / 3) * 3;

        double S = 0, N = 0;   // Potential sites
        double Sd = 0, Nd = 0; // Observed differences

        for (int i = 0; i < len; i += 3) {
            String codon1 = seq1.substring(i, i + 3);
            String codon2 = seq2.substring(i, i + 3);

            if (codon1.contains("N") || codon2.contains("N") || 
                codon1.contains("-") || codon2.contains("-")) continue;

            // 1. Calculate potential S and N for codon1 and codon2 (average)
            double[] sn1 = countPotentialSites(codon1);
            double[] sn2 = countPotentialSites(codon2);
            S += (sn1[0] + sn2[0]) / 2.0;
            N += (sn1[1] + sn2[1]) / 2.0;

            // 2. Calculate differences
            if (!codon1.equals(codon2)) {
                double[] diffs = countDifferences(codon1, codon2);
                Sd += diffs[0];
                Nd += diffs[1];
            }
        }

        double pS = (S > 0) ? Sd / S : 0;
        double pN = (N > 0) ? Nd / N : 0;

        // Jukes-Cantor correction
        double ks = (pS < 0.75) ? -0.75 * Math.log(1 - (4.0/3.0) * pS) : pS;
        double ka = (pN < 0.75) ? -0.75 * Math.log(1 - (4.0/3.0) * pN) : pN;

        return new Result(ka, ks);
    }

    private double[] countPotentialSites(String codon) {
        double s = 0, n = 0;
        String aa = GENETIC_CODE.getOrDefault(codon, "X");
        
        for (int i = 0; i < 3; i++) {
            int synCount = 0;
            char original = codon.charAt(i);
            for (char base : new char[]{'A', 'C', 'G', 'T'}) {
                if (base == original) continue;
                StringBuilder sb = new StringBuilder(codon);
                sb.setCharAt(i, base);
                String mutAA = GENETIC_CODE.getOrDefault(sb.toString(), "X");
                if (mutAA.equals(aa)) synCount++;
            }
            s += synCount / 3.0;
            n += (3 - synCount) / 3.0;
        }
        return new double[]{s, n};
    }

    private double[] countDifferences(String c1, String c2) {
        int diffs = 0;
        for (int i = 0; i < 3; i++) if (c1.charAt(i) != c2.charAt(i)) diffs++;
        
        if (diffs == 1) {
            String aa1 = GENETIC_CODE.getOrDefault(c1, "X");
            String aa2 = GENETIC_CODE.getOrDefault(c2, "X");
            return aa1.equals(aa2) ? new double[]{1, 0} : new double[]{0, 1};
        } else {
            // For 2 or 3 diffs, compute all paths and average (simplified for this tool)
            // Here we use a simpler approach: proportional to potential sites
            double[] sn1 = countPotentialSites(c1);
            double total = sn1[0] + sn1[1];
            return new double[]{diffs * (sn1[0]/total), diffs * (sn1[1]/total)};
        }
    }

    /**
     * Basic Needleman-Wunsch for global alignment of CDS.
     */
    public String[] align(String s1, String s2) {
        int match = 2, mismatch = -1, gap = -2;
        int n = s1.length(), m = s2.length();
        int[][] dp = new int[n + 1][m + 1];

        for (int i = 0; i <= n; i++) dp[i][0] = i * gap;
        for (int j = 0; j <= m; j++) dp[0][j] = j * gap;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int score = (s1.charAt(i - 1) == s2.charAt(j - 1)) ? match : mismatch;
                dp[i][j] = Math.max(dp[i - 1][j - 1] + score, Math.max(dp[i - 1][j] + gap, dp[i][j - 1] + gap));
            }
        }

        StringBuilder a1 = new StringBuilder(), a2 = new StringBuilder();
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + ((s1.charAt(i - 1) == s2.charAt(j - 1)) ? match : mismatch)) {
                a1.append(s1.charAt(i - 1)); a2.append(s2.charAt(j - 1)); i--; j--;
            } else if (i > 0 && dp[i][j] == dp[i - 1][j] + gap) {
                a1.append(s1.charAt(i - 1)); a2.append('-'); i--;
            } else {
                a1.append('-'); a2.append(s2.charAt(j - 1)); j--;
            }
        }
        return new String[]{a1.reverse().toString(), a2.reverse().toString()};
    }
}
