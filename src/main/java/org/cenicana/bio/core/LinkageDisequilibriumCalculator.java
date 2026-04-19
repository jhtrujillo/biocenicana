package org.cenicana.bio.core;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;
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

    public void computeLD(String inputVcf, String outputTsv) throws IOException {
        String[] sampleNames = VcfFastReader.getSampleIds(inputVcf);
        int numSamples = sampleNames.length;

        System.out.println("[LD] Computing Linkage Disequilibrium with allele dosages...");
        System.out.println("[LD] Max distance: " + maxDistanceBp + " bp. Min MAF: " + minMaf);

        try (BufferedReader br = new BufferedReader(new FileReader(inputVcf));
             PrintWriter pw = new PrintWriter(new FileWriter(outputTsv))) {

            // Header for TSV
            pw.println("CHR\tPOS1\tID1\tPOS2\tID2\tDISTANCE\tR2");

            // Bins for LD Decay Curve
            int numBins = (maxDistanceBp / binSizeBp) + 1;
            double[] sumR2 = new double[numBins];
            long[] countR2 = new long[numBins];

            LinkedList<MarkerDosage> slidingWindow = new LinkedList<>();

            String line;
            int totalProcessed = 0;
            int totalPairs = 0;

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

                for (int i = 0; i < numSamples; i++) {
                    String[] gData = cols[9 + i].split(":");
                    double countRef = 0;
                    double countAlt = 0;
                    boolean foundCounts = false;

                    try {
                        if (bsdpIdx != -1 && gData.length > bsdpIdx) {
                            String[] bsdp = gData[bsdpIdx].split(",");
                            if (bsdp.length >= 2 && !bsdp[0].equals(".") && !bsdp[1].equals(".")) {
                                countAlt = Double.parseDouble(bsdp[0]);
                                countRef = Double.parseDouble(bsdp[1]);
                                foundCounts = true;
                            }
                        } else if (adIdx != -1 && gData.length > adIdx) {
                            String[] ads = gData[adIdx].split(",");
                            if (ads.length >= 2 && !ads[0].equals(".") && !ads[1].equals(".")) {
                                countRef = Double.parseDouble(ads[0]);
                                countAlt = Double.parseDouble(ads[1]);
                                foundCounts = true;
                            }
                        } else if (roIdx != -1 && aoIdx != -1 && gData.length > Math.max(roIdx, aoIdx)) {
                            if (!gData[roIdx].equals(".") && !gData[aoIdx].equals(".")) {
                                countRef = Double.parseDouble(gData[roIdx]);
                                String aoStr = gData[aoIdx].contains(",") ? gData[aoIdx].split(",")[0] : gData[aoIdx];
                                countAlt = Double.parseDouble(aoStr);
                                foundCounts = true;
                            }
                        }
                    } catch (NumberFormatException ignored) {}

                    if (foundCounts && (countRef + countAlt) >= 5) { // Min 5 reads depth to trust dosage
                        // Dosage = proportion of alternate allele
                        double dosage = countAlt / (countRef + countAlt);
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

                // Calculate pairwise LD with all markers currently in the window
                for (MarkerDosage pastMarker : slidingWindow) {
                    double r2 = calculatePearsonR2(marker, pastMarker, numSamples);
                    if (r2 >= 0.0) { // Valid correlation calculated
                        int distance = pos - pastMarker.pos;
                        pw.printf("%s\t%d\t%s\t%d\t%s\t%d\t%.6f%n",
                                chrom, pastMarker.pos, pastMarker.id, pos, id, distance, r2);
                        
                        // Accumulate for LD Decay
                        if (distance <= maxDistanceBp) {
                            int binIdx = distance / binSizeBp;
                            sumR2[binIdx] += r2;
                            countR2[binIdx]++;
                        }
                        
                        totalPairs++;
                    }
                }

                slidingWindow.addLast(marker);
                totalProcessed++;
            }

            System.out.println("[LD] Complete. Processed markers: " + totalProcessed);
            System.out.println("[LD] Pairwise LD combinations calculated: " + totalPairs);

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
