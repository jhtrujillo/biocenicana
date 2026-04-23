package org.cenicana.bio.core;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
        public double[] allDosages; // Individual sample dosages (0.0 - 1.0 frequency)
        public int[] refDepths;     // Reference allele depths
        public int[] altDepths;     // Alternative allele depths
        public float minVal, maxVal;
    }

    public static class SampleCoord {
        public String name;
        public double x, y;
    }

    private int numBins = 50; 
    private int maxClusters = 11; // Standard for ploidy 10
    private String[] sampleNames; 

    public List<SnpResult> analyzeVcf(String vcfFile, int ploidy) {
        AlleleDosageCalculator dosageCalc = new AlleleDosageCalculator();
        dosageCalc.maxSnps = 1000; 
        
        List<AlleleDosageCalculator.DosageResult> vcfResults = dosageCalc.calculate(vcfFile, ploidy);
        
        List<SnpResult> snpResults = new ArrayList<>();
        for (AlleleDosageCalculator.DosageResult dr : vcfResults) {
            SnpResult sr = new SnpResult();
            sr.id = dr.snpId;
            sr.allDosages = dr.dosages;
            sr.refDepths = dr.refDepths;
            sr.altDepths = dr.altDepths;
            
            float[] validData = new float[dr.dosages.length];
            int count = 0;
            for (double d : dr.dosages) {
                if (d >= 0) validData[count++] = (float)d;
            }
            float[] data = Arrays.copyOf(validData, count);
            
            sr.histogramBins = calculateHistogram(data);
            performClustering(sr, data, ploidy);
            snpResults.add(sr);
        }
        return snpResults;
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
                res.allDosages = new double[cols.length - 2];
                int count = 0;

                for (int i = 2; i < cols.length; i++) {
                    float val = Float.parseFloat(cols[i]);
                    res.allDosages[i-2] = val;
                    if (val >= 0) {
                        dosages[count++] = val;
                    }
                }

                if (count == 0) continue;
                float[] validDosages = Arrays.copyOf(dosages, count);
                res.histogramBins = calculateHistogram(validDosages);
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
        int k = Math.min(dosages.length, ploidy + 1);
        float[] centroids = new float[ploidy + 1];
        for (int i = 0; i <= ploidy; i++) centroids[i] = (1.0f / ploidy) * i;

        int[] assignments = new int[dosages.length];
        boolean changed = true;
        int maxIter = 20;

        for (int iter = 0; iter < maxIter && changed; iter++) {
            changed = false;
            for (int i = 0; i < dosages.length; i++) {
                float bestDist = Float.MAX_VALUE;
                int bestK = 0;
                for (int j = 0; j <= ploidy; j++) {
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
            float[] sums = new float[ploidy + 1];
            int[] counts = new int[ploidy + 1];
            for (int i = 0; i < dosages.length; i++) {
                sums[assignments[i]] += dosages[i];
                counts[assignments[i]]++;
            }
            for (int j = 0; j <= ploidy; j++) {
                if (counts[j] > 0) centroids[j] = sums[j] / counts[j];
            }
        }

        res.centroids = centroids;
        res.clusterCounts = new int[ploidy + 1];
        for (int i = 0; i < dosages.length; i++) {
            res.clusterCounts[assignments[i]]++;
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
                if(hCols[i].equalsIgnoreCase("PC1")) pc1Idx = i;
                if(hCols[i].equalsIgnoreCase("PC2")) pc2Idx = i;
            }
            if(pc1Idx == -1 || pc2Idx == -1) return coords;

            String line;
            while ((line = br.readLine()) != null) {
                String[] cols = line.split(",");
                if(cols.length <= Math.max(pc1Idx, pc2Idx)) continue;
                SampleCoord sc = new SampleCoord();
                sc.name = cols[0].replace("\"", "");
                sc.x = Double.parseDouble(cols[pc1Idx]);
                sc.y = Double.parseDouble(cols[pc2Idx]);
                coords.add(sc);
            }
        }
        return coords;
    }
}
