package org.cenicana.bio.core;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Analyzes SNP dosage distributions by performing 1D clustering and histogram calculation.
 * This is used by the SNP Quality Explorer to visualize how dosages group across the population.
 */
public class SnpClusterAnalyzer {

    public static class SnpResult {
        public String id;
        public String chr;
        public int pos;
        public float[] centroids;
        public int[] clusterCounts;
        public int[] histogramBins;
        public float minVal, maxVal;
    }

    private int numBins = 25;
    private int maxClusters = 11; // Standard for ploidy 10

    public List<SnpResult> analyzeMatrix(String matrixFile, int ploidy) throws IOException {
        List<SnpResult> results = new ArrayList<>();
        this.maxClusters = ploidy + 1;

        try (BufferedReader br = new BufferedReader(new FileReader(matrixFile))) {
            String header = br.readLine();
            if (header == null) return results;

            String line;
            while ((line = br.readLine()) != null) {
                String[] cols = line.split("\t");
                if (cols.length < 3) continue;

                SnpResult res = new SnpResult();
                res.chr = cols[0];
                res.pos = Integer.parseInt(cols[1]);
                res.id = res.chr + "_" + res.pos;

                float[] dosages = new float[cols.length - 2];
                int count = 0;
                res.minVal = 1.0f;
                res.maxVal = 0.0f;

                for (int i = 2; i < cols.length; i++) {
                    float val = Float.parseFloat(cols[i]);
                    if (val >= 0) {
                        dosages[count++] = val;
                        if (val < res.minVal) res.minVal = val;
                        if (val > res.maxVal) res.maxVal = val;
                    }
                }

                if (count == 0) continue;
                float[] validDosages = Arrays.copyOf(dosages, count);

                // 1. Calculate Histogram
                res.histogramBins = calculateHistogram(validDosages);

                // 2. Perform 1D Clustering (K-Means)
                performClustering(res, validDosages, ploidy);

                results.add(res);
            }
        }
        return results;
    }

    private int[] calculateHistogram(float[] dosages) {
        int[] bins = new int[numBins];
        for (float d : dosages) {
            int idx = (int) (d * (numBins - 1));
            if (idx >= 0 && idx < numBins) bins[idx]++;
        }
        return bins;
    }

    private void performClustering(SnpResult res, float[] dosages, int ploidy) {
        int k = Math.min(dosages.length, maxClusters);
        float[] centroids = new float[k];
        
        // Initial centroids: theoretical ploidy levels
        for (int i = 0; i < k; i++) {
            centroids[i] = (1.0f / ploidy) * i;
        }

        int[] assignments = new int[dosages.length];
        boolean changed = true;
        int maxIter = 20;

        for (int iter = 0; iter < maxIter && changed; iter++) {
            changed = false;
            // E-step: Assign to nearest
            for (int i = 0; i < dosages.length; i++) {
                float bestDist = Float.MAX_VALUE;
                int bestK = 0;
                for (int j = 0; j < k; j++) {
                    float dist = Math.abs(dosages[i] - centroids[j]);
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestK = j;
                    }
                }
                if (assignments[i] != bestK) {
                    assignments[i] = bestK;
                    changed = true;
                }
            }

            // M-step: Update centroids
            float[] sums = new float[k];
            int[] counts = new int[k];
            for (int i = 0; i < dosages.length; i++) {
                sums[assignments[i]] += dosages[i];
                counts[assignments[i]]++;
            }
            for (int j = 0; j < k; j++) {
                if (counts[j] > 0) centroids[j] = sums[j] / counts[j];
            }
        }

        // Filter out empty clusters and save
        List<Float> finalCentroids = new ArrayList<>();
        List<Integer> finalCounts = new ArrayList<>();
        for (int j = 0; j < k; j++) {
            int count = 0;
            float sum = 0;
            for(int i=0; i<dosages.length; i++) {
                if(assignments[i] == j) {
                    count++;
                    sum += dosages[i];
                }
            }
            if (count > 0) {
                finalCentroids.add(sum / count);
                finalCounts.add(count);
            }
        }

        res.centroids = new float[finalCentroids.size()];
        res.clusterCounts = new int[finalCounts.size()];
        for (int i = 0; i < finalCentroids.size(); i++) {
            res.centroids[i] = finalCentroids.get(i);
            res.clusterCounts[i] = finalCounts.get(i);
        }
    }
}
