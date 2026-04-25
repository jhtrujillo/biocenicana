package org.cenicana.bio.core;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.cenicana.bio.io.VcfFastReader;

/**
 * Calculates Linkage Disequilibrium (LD) r-squared values between markers 
 * using allele dosages (continuous frequencies) to support polyploid crops like Sugarcane.
 * Automatically handles VCFs from NGSEP, FreeBayes, and GATK.
 */
public class LinkageDisequilibriumCalculator {

    private int ploidy = 2;
    private double minMaf = 0.05;
    private int maxDistanceBp = 50000;
    private boolean generateHtml = false;
    private int binSizeBp = 1000; // Default bin size for decay curve

    public void setPloidy(int ploidy) { this.ploidy = ploidy; }
    public void setMinMaf(double minMaf) { this.minMaf = minMaf; }
    public void setMaxDistanceBp(int maxDistanceBp) { this.maxDistanceBp = maxDistanceBp; }
    public void setGenerateHtml(boolean generateHtml) { this.generateHtml = generateHtml; }
    public void setBinSizeBp(int binSizeBp) { this.binSizeBp = binSizeBp; }

    private static class MarkerDosage {
        String chrom;
        int pos;
        String id;
        double[] dosages;
        boolean[] missing;
        double maf;

        MarkerDosage(String chrom, int pos, String id, int numSamples) {
            this.chrom = chrom;
            this.pos = pos;
            this.id = id;
            this.dosages = new double[numSamples];
            this.missing = new boolean[numSamples];
        }
    }

    public void computeLD(String inputVcf, String outputTsv, int threads) throws IOException {
        String[] sampleNames = VcfFastReader.getSampleIds(inputVcf);
        int numSamples = sampleNames.length;
        
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        final AtomicInteger totalPairs = new AtomicInteger(0);
        final AtomicInteger totalProcessed = new AtomicInteger(0);

        System.out.println("[LD] Computing Linkage Disequilibrium with allele dosages...");
        System.out.println("[LD] Max distance: " + maxDistanceBp + " bp. Min MAF: " + minMaf);

        try (BufferedReader br = new BufferedReader(new FileReader(inputVcf));
             PrintWriter pw = new PrintWriter(new FileWriter(outputTsv))) {

            // Header for TSV
            pw.println("CHR\tPOS1\tID1\tPOS2\tID2\tDISTANCE\tR2");

            // Bins for LD Decay Curve
            int numBins = (maxDistanceBp / binSizeBp) + 1;
            final double[] sumR2 = new double[numBins];
            final long[] countR2 = new long[numBins];
            final Object matrixLock = new Object();

            LinkedList<MarkerDosage> slidingWindow = new LinkedList<>();

            String line;

            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] cols = line.split("\t");
                if (cols.length < 9 + numSamples) continue;

                String chrom = cols[0];
                int pos = Integer.parseInt(cols[1]);
                String id = cols[2];
                String format = cols[8];

                int adIdx = -1, roIdx = -1, aoIdx = -1, bsdpIdx = -1;
                String[] formatTokens = format.split(":");
                for (int i = 0; i < formatTokens.length; i++) {
                    switch (formatTokens[i]) {
                        case "AD":   adIdx   = i; break;
                        case "RO":   roIdx   = i; break;
                        case "AO":   aoIdx   = i; break;
                        case "BSDP": bsdpIdx = i; break;
                    }
                }

                MarkerDosage marker = new MarkerDosage(chrom, pos, id, numSamples);
                int genotypedCount = 0;
                double sumDosage = 0;
                String refBase = cols[3];
                String altBase = cols[4];

                for (int i = 0; i < numSamples; i++) {
                    String[] gData = cols[9 + i].split(":");

                    double[] counts = AlleleDosageCalculator.getRefAltCounts(gData, adIdx, bsdpIdx, roIdx, aoIdx, refBase, altBase);

                    if (counts != null && (counts[0] + counts[1]) >= 5) {
                        double dosage = counts[1] / (counts[0] + counts[1]);
                        marker.dosages[i] = dosage;
                        marker.missing[i] = false;
                        genotypedCount++;
                        sumDosage += dosage;
                    } else {
                        marker.missing[i] = true;
                    }
                }

                // Filter by MAF
                if (genotypedCount == 0) continue;
                double meanDosage = sumDosage / genotypedCount;
                double maf = Math.min(meanDosage, 1.0 - meanDosage);
                marker.maf = maf;

                if (maf < minMaf) continue;

                // Remove markers from sliding window that are on a different chromosome or too far away
                while (!slidingWindow.isEmpty()) {
                    MarkerDosage first = slidingWindow.peekFirst();
                    if (!first.chrom.equals(chrom) || (pos - first.pos) > maxDistanceBp) {
                        slidingWindow.pollFirst();
                    } else {
                        break; // Sorted by pos, so if first is within range, rest are too
                    }
                }

                // Calculate pairwise LD with all markers currently in the window in parallel
                if (!slidingWindow.isEmpty()) {
                    List<MarkerDosage> currentWindow = new ArrayList<>(slidingWindow);
                    final MarkerDosage currentMarker = marker;
                    
                    executor.submit(() -> {
                        for (MarkerDosage pastMarker : currentWindow) {
                            double r2 = calculatePearsonR2(currentMarker, pastMarker, numSamples);
                            if (r2 >= 0.0) {
                                int distance = currentMarker.pos - pastMarker.pos;
                                synchronized (pw) {
                                    pw.printf("%s\t%d\t%s\t%d\t%s\t%d\t%.6f%n",
                                            currentMarker.chrom, pastMarker.pos, pastMarker.id, currentMarker.pos, currentMarker.id, distance, r2);
                                }
                                
                                synchronized (matrixLock) {
                                    int binIdx = distance / binSizeBp;
                                    if (binIdx < numBins) {
                                        sumR2[binIdx] += r2;
                                        countR2[binIdx]++;
                                    }
                                }
                                totalPairs.incrementAndGet();
                            }
                        }
                    });
                }

                slidingWindow.addLast(marker);
                totalProcessed.incrementAndGet();
            }

            executor.shutdown();
            try { executor.awaitTermination(1, TimeUnit.DAYS); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }

            System.out.println("[LD] Complete. Processed markers: " + totalProcessed.get());
            System.out.println("[LD] Pairwise LD combinations calculated: " + totalPairs.get());

            // Calculate LD Half-Decay Distance
            double maxR2 = 0.0;
            int halfDecayDistanceBp = -1;
            
            for (int i = 0; i < sumR2.length; i++) {
                if (countR2[i] > 0) {
                    double avg = sumR2[i] / countR2[i];
                    if (avg > maxR2) maxR2 = avg;
                }
            }
            
            double thresholdR2 = maxR2 / 2.0;
            for (int i = 0; i < sumR2.length; i++) {
                if (countR2[i] > 0) {
                    double avg = sumR2[i] / countR2[i];
                    if (avg <= thresholdR2) {
                        halfDecayDistanceBp = i * binSizeBp;
                        break;
                    }
                }
            }

            System.out.println("[LD] Maximum Average R2 (Baseline): " + String.format(java.util.Locale.US, "%.4f", maxR2));
            if (halfDecayDistanceBp != -1) {
                System.out.println("[LD] Half-Decay Distance: ~" + halfDecayDistanceBp + " bp  (R2 drops below " + String.format(java.util.Locale.US, "%.4f", thresholdR2) + ")");
            } else {
                System.out.println("[LD] Half-Decay Distance: Not reached within " + maxDistanceBp + " bp. Consider increasing window size.");
            }

            if (generateHtml) {
                System.out.println("[LD] Generating Interactive LD Decay Dashboard...");
                String htmlFile = outputTsv.endsWith(".tsv") ? outputTsv.replace(".tsv", "_decay.html") : outputTsv + "_decay.html";
                org.cenicana.bio.io.HtmlDashboardGenerator.generateLdDecayDashboard(htmlFile, sumR2, countR2, binSizeBp, halfDecayDistanceBp, thresholdR2);
                System.out.println("[LD] HTML Dashboard generated at: " + htmlFile);
            }
        }
    }

    private double calculatePearsonR2(MarkerDosage m1, MarkerDosage m2, int numSamples) {
        double sum1 = 0, sum2 = 0, sum1Sq = 0, sum2Sq = 0, pSum = 0;
        int n = 0;

        for (int i = 0; i < numSamples; i++) {
            if (!m1.missing[i] && !m2.missing[i]) {
                double x = m1.dosages[i];
                double y = m2.dosages[i];
                sum1 += x;
                sum2 += y;
                sum1Sq += x * x;
                sum2Sq += y * y;
                pSum += x * y;
                n++;
            }
        }

        if (n < 3) return -1.0; // Need at least 3 common valid samples to compute correlation

        double num = pSum - (sum1 * sum2 / n);
        double den = Math.sqrt((sum1Sq - (sum1 * sum1) / n) * (sum2Sq - (sum2 * sum2) / n));

        if (den == 0) return 0.0;
        double r = num / den;
        return r * r;
    }
}
