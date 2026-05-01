package org.cenicana.bio.core;

import java.io.*;
import java.util.*;

/**
 * High-performance VCF merger using Multi-way Merge algorithm.
 * Streams data from multiple files simultaneously to join them by coordinates.
 * Samples missing in a source file are filled with "./.".
 */
public class VcfMerger {

    private List<String> inputFiles;
    private File outputFile;

    public VcfMerger(List<String> inputFiles, String outputPath) {
        this.inputFiles = inputFiles;
        this.outputFile = new File(outputPath);
    }

    public void merge() throws IOException {
        List<BufferedReader> readers = new ArrayList<>();
        List<String> allSamples = new ArrayList<>();
        Map<Integer, List<String>> fileToSamples = new HashMap<>();

        // 1. Initial Scan for Samples and Headers
        System.out.println("[Vcf-Merge] Collecting sample headers from " + inputFiles.size() + " files...");
        for (int i = 0; i < inputFiles.size(); i++) {
            try (BufferedReader br = new BufferedReader(new FileReader(inputFiles.get(i)))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#CHROM")) {
                        String[] cols = line.split("\t");
                        List<String> samples = new ArrayList<>();
                        for (int j = 9; j < cols.length; j++) {
                            samples.add(cols[j]);
                            if (!allSamples.contains(cols[j])) {
                                allSamples.add(cols[j]);
                            }
                        }
                        fileToSamples.put(i, samples);
                        break;
                    }
                }
            }
        }
        Collections.sort(allSamples);
        System.out.println("[Vcf-Merge] Total unique samples found: " + allSamples.size());

        // 2. Setup Multi-way stream
        PriorityQueue<VcfLine> pq = new PriorityQueue<>();
        for (int i = 0; i < inputFiles.size(); i++) {
            BufferedReader br = new BufferedReader(new FileReader(inputFiles.get(i)));
            readers.add(br);
            VcfLine firstLine = readNextGenomicLine(br, i);
            if (firstLine != null) pq.add(firstLine);
        }

        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)))) {
            // Write Header
            writer.println("##fileformat=VCFv4.2");
            writer.println("##source=BioJavaVcfMerger");
            writer.println("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples with Data\">");
            writer.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            writer.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for (String sample : allSamples) writer.print("\t" + sample);
            writer.println();

            // 3. Process Streams
            while (!pq.isEmpty()) {
                List<VcfLine> currentSites = new ArrayList<>();
                VcfLine top = pq.poll();
                currentSites.add(top);

                // Pull all lines matching the same coordinate (Chr + Pos)
                while (!pq.isEmpty() && pq.peek().compareTo(top) == 0) {
                    currentSites.add(pq.poll());
                }

                // Merge genotypes for this site
                VcfLine lead = currentSites.get(0);
                StringBuilder mergedLine = new StringBuilder();
                mergedLine.append(lead.chr).append("\t")
                          .append(lead.pos).append("\t")
                          .append(lead.id).append("\t")
                          .append(lead.ref).append("\t")
                          .append(lead.alt).append("\t.\t.\t.\tGT");

                // Fill genotypes for all samples
                for (String sample : allSamples) {
                    String genotype = "./.";
                    for (VcfLine site : currentSites) {
                        List<String> sourceSamples = fileToSamples.get(site.fileIndex);
                        int sampleIdx = sourceSamples.indexOf(sample);
                        if (sampleIdx != -1) {
                            genotype = site.genotypes[sampleIdx];
                            break; 
                        }
                    }
                    mergedLine.append("\t").append(genotype);
                }
                writer.println(mergedLine.toString());

                // Refill PQ from the source files we just used
                for (VcfLine site : currentSites) {
                    VcfLine next = readNextGenomicLine(readers.get(site.fileIndex), site.fileIndex);
                    if (next != null) pq.add(next);
                }
            }
        } finally {
            for (BufferedReader br : readers) br.close();
        }
    }

    private VcfLine readNextGenomicLine(BufferedReader br, int fileIndex) throws IOException {
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) continue;
            return new VcfLine(line, fileIndex);
        }
        return null;
    }

    /**
     * Helper class to represent a VCF line in the PriorityQueue
     */
    private static class VcfLine implements Comparable<VcfLine> {
        String chr;
        int pos;
        String id, ref, alt;
        String[] genotypes;
        int fileIndex;

        public VcfLine(String line, int fileIndex) {
            String[] cols = line.split("\t");
            this.chr = cols[0];
            this.pos = Integer.parseInt(cols[1]);
            this.id = cols[2];
            this.ref = cols[3];
            this.alt = cols[4];
            this.fileIndex = fileIndex;
            this.genotypes = new String[cols.length - 9];
            System.arraycopy(cols, 9, this.genotypes, 0, cols.length - 9);
        }

        @Override
        public int compareTo(VcfLine other) {
            int c = this.chr.compareTo(other.chr);
            if (c != 0) return c;
            return Integer.compare(this.pos, other.pos);
        }
    }
}
