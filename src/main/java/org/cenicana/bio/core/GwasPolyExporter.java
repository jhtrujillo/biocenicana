package org.cenicana.bio.core;

import org.cenicana.bio.io.VcfFastReader;
import java.io.*;
import java.util.*;

/**
 * Exporter for GWASpoly format.
 * Converts allele dosages to ACGT strings (e.g., AAAAAAAAAA or AACCCCCCCC for ploidy 10).
 */
public class GwasPolyExporter {

    public void exportToACGT(String vcfFile, String outputFile, int ploidy) throws IOException {
        String[] sampleNames = VcfFastReader.getSampleIds(vcfFile);
        int numSamples = sampleNames.length;

        try (BufferedReader br = new BufferedReader(new FileReader(vcfFile));
             PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {

            // Header: Marker, Chrom, Pos, [Samples...]
            pw.print("Marker,Chrom,Pos");
            for (String sample : sampleNames) {
                pw.print("," + sample);
            }
            pw.println();

            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] cols = line.split("\t");
                if (cols.length < 9 + numSamples) continue;

                String marker = cols[2].equals(".") ? cols[0] + "_" + cols[1] : cols[2];
                String chrom  = cols[0];
                String pos    = cols[1];
                String ref    = cols[3];
                String alt    = cols[4];

                // We only support biallelic SNPs for this simple exporter
                if (ref.length() > 1 || alt.length() > 1 || alt.contains(",")) continue;

                pw.print(marker + "," + chrom + "," + pos);

                // Identify depth indices
                String[] format = cols[8].split(":");
                int adIdx = -1, bsdpIdx = -1, roIdx = -1, aoIdx = -1;
                for (int i = 0; i < format.length; i++) {
                    switch (format[i]) {
                        case "AD":   adIdx   = i; break;
                        case "BSDP": bsdpIdx = i; break;
                        case "RO":   roIdx   = i; break;
                        case "AO":   aoIdx   = i; break;
                    }
                }

                for (int i = 0; i < numSamples; i++) {
                    String sampleData = cols[9 + i];
                    String[] fields = sampleData.split(":");
                    
                    String acgtString = "";
                    double dosage = -1; // -1 means missing

                    // 1. Try to get dosage from read counts using centralized parser
                    double[] counts = AlleleDosageCalculator.getRefAltCounts(fields, adIdx, bsdpIdx, roIdx, aoIdx, ref, alt);
                    if (counts != null && counts[0] + counts[1] > 0) {
                        dosage = counts[0] / (counts[0] + counts[1]); // ref proportion
                    }

                    // 2. Fallback to GT if AD is missing
                    if (dosage == -1) {
                        String gt = fields[0];
                        if (!gt.equals(".") && !gt.startsWith("./")) {
                            String[] alleles = gt.split("[/|]");
                            int rCount = 0;
                            int total = alleles.length;
                            for (String a : alleles) {
                                if (a.equals("0")) rCount++;
                            }
                            if (total > 0) {
                                dosage = (double) rCount / total;
                            }
                        }
                    }

                    // 3. Generate string
                    if (dosage == -1) {
                        // Missing: ploidy times 'N'
                        char[] nArr = new char[ploidy];
                        Arrays.fill(nArr, 'N');
                        acgtString = new String(nArr);
                    } else {
                        int nRef = (int) Math.round(dosage * ploidy);
                        int nAlt = ploidy - nRef;
                        
                        StringBuilder sb = new StringBuilder();
                        for (int j = 0; j < nRef; j++) sb.append(ref);
                        for (int j = 0; j < nAlt; j++) sb.append(alt);
                        acgtString = sb.toString();
                    }

                    pw.print("," + acgtString);
                }
                pw.println();
            }
        }
        System.out.println("[GwasPolyExporter] Export complete: " + outputFile);
    }
}
