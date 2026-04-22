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
        public float[] allDosages; // New: Full population dosages
        public float minVal, maxVal;
    }

    public static class SampleCoord {
        public String name;
        public double x, y;
    }

    private int numBins = 25;
    private int maxClusters = 11; // Standard for ploidy 10
    private String[] sampleNames; // To keep track of sample order
    private List<SampleCoord> pcaCoords;

    public void setPcaCoords(List<SampleCoord> coords) {
        this.pcaCoords = coords;
    }

    public List<SnpResult> analyzeMatrix(String matrixFile, int ploidy) throws IOException {
        List<SnpResult> results = new ArrayList<>();
        this.maxClusters = ploidy + 1;

        try (BufferedReader br = new BufferedReader(new FileReader(matrixFile))) {
            String header = br.readLine();
            if (header == null) return results;
            String[] headerCols = header.split("\t");
            if (headerCols.length > 2) {
                this.sampleNames = Arrays.copyOfRange(headerCols, 2, headerCols.length);
            }

            String line;
            while ((line = br.readLine()) != null) {
                String[] cols = line.split("\t");
                if (cols.length < 3) continue;

                SnpResult res = new SnpResult();
                res.chr = cols[0];
                res.pos = Integer.parseInt(cols[1]);
                res.id = res.chr + "_" + res.pos;

                float[] dosages = new float[cols.length - 2];
                res.allDosages = new float[cols.length - 2]; // Store all for PCA
                int count = 0;
                res.minVal = 1.0f;
                res.maxVal = 0.0f;

                for (int i = 2; i < cols.length; i++) {
                    float val = Float.parseFloat(cols[i]);
                    res.allDosages[i-2] = val; // Even if -1
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

    public List<SampleCoord> loadPcaCsv(String pcaCsv) throws IOException {
        List<SampleCoord> coords = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(pcaCsv))) {
            String header = br.readLine();
            if (header == null) return coords;
            String[] hCols = header.split(",");
            int pc1Idx = -1, pc2Idx = -1;
            for(int i=0; i<hCols.length; i++) {
                if(hCols[i].equals("PC1")) pc1Idx = i;
                if(hCols[i].equals("PC2")) pc2Idx = i;
            }
            if(pc1Idx == -1 || pc2Idx == -1) return coords;

            String line;
            while ((line = br.readLine()) != null) {
                String[] cols = line.split(",");
                if(cols.length <= Math.max(pc1Idx, pc2Idx)) continue;
                SampleCoord sc = new SampleCoord();
                sc.name = cols[0];
                sc.x = Double.parseDouble(cols[pc1Idx]);
                sc.y = Double.parseDouble(cols[pc2Idx]);
                coords.add(sc);
            }
        }
        return coords;
    }

    public String[] getSampleNames() {
        return sampleNames;
    }
}
