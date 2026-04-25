package org.cenicana.bio.utils;

/**
 * Utility class for Hardy-Weinberg Equilibrium (HWE) statistical tests.
 */
public class HweUtils {

    /**
     * Calculates the Fisher Exact Test p-value for HWE in a biallelic site (diploid).
     * Based on the algorithm by Wigginton, Cutler, and Abecasis (2005).
     *
     * @param obsHom1 Number of homozygous reference genotypes
     * @param obsHet  Number of heterozygous genotypes
     * @param obsHom2 Number of homozygous alternative genotypes
     * @return The p-value (0 to 1) for the null hypothesis of HWE.
     */
    public static double calculateHweFisherExactTest(int obsHom1, int obsHet, int obsHom2) {
        int n = obsHom1 + obsHet + obsHom2;
        int nA = 2 * obsHom1 + obsHet;
        int nB = 2 * obsHom2 + obsHet;

        int minorAlleleCount = Math.min(nA, nB);
        int majorAlleleCount = Math.max(nA, nB);

        int minHet = minorAlleleCount % 2;
        int maxHet = minorAlleleCount;

        double[] probs = new double[maxHet + 1];

        // Using double for logFacs to avoid overflow and maintain precision
        double[] logFacs = new double[n * 2 + 1];
        for (int i = 1; i <= n * 2; i++) {
            logFacs[i] = logFacs[i - 1] + Math.log(i);
        }

        double maxLogProb = -Double.MAX_VALUE;
        for (int h = minHet; h <= maxHet; h += 2) {
            int h1 = (minorAlleleCount - h) / 2;
            int h2 = n - h - h1;
            double logP = logFacs[n] + logFacs[minorAlleleCount] + logFacs[majorAlleleCount]
                - logFacs[2 * n] - logFacs[h1] - logFacs[h] - logFacs[h2] + h * Math.log(2.0);
            probs[h] = logP;
            if (logP > maxLogProb) {
                maxLogProb = logP;
            }
        }

        double sum = 0.0;
        for (int h = minHet; h <= maxHet; h += 2) {
            probs[h] = Math.exp(probs[h] - maxLogProb);
            sum += probs[h];
        }

        double pValue = 0.0;
        double obsProb = probs[obsHet];
        for (int h = minHet; h <= maxHet; h += 2) {
            if (probs[h] <= obsProb + 1e-9) {
                pValue += probs[h] / sum;
            }
        }
        return Math.min(1.0, pValue);
    }
}
