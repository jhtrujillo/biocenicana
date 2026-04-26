package org.cenicana.bio.core;

import org.cenicana.bio.io.FastaReader;
import org.cenicana.bio.model.Gene;
import org.cenicana.bio.model.SyntenicBlock;
import org.cenicana.bio.model.SyntenicPair;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Phylogeny Module for BioCenicana.
 * Calculates evolutionary distances and prepares data for tree construction.
 */
public class PhylogenyAnalyzer {

    public void runPhylogeny(List<SyntenicBlock> blocks, Map<String, Gene> genes1, Map<String, Gene> genes2,
                             String fasta1, String fasta2, String outputPrefix) throws IOException {
        
        System.out.println("[Phylogeny] Loading sequences for orthologs...");
        FastaReader reader = new FastaReader();
        Map<String, String> seqs1 = reader.read(fasta1);
        Map<String, String> seqs2 = reader.read(fasta2);

        System.out.println("[Phylogeny] Calculating evolutionary distances...");
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputPrefix + "_ortholog_distances.tsv"))) {
            pw.println("Gene1\tGene2\tIdentity\tJukesCantorDistance");

            for (SyntenicBlock block : blocks) {
                for (SyntenicPair pair : block.getPairs()) {
                    String s1 = seqs1.get(pair.getGeneId1());
                    String s2 = seqs2.get(pair.getGeneId2());

                    if (s1 != null && s2 != null) {
                        double dist = calculateDistance(s1, s2);
                        double identity = calculateIdentity(s1, s2);
                        pw.printf(Locale.US, "%s\t%s\t%.4f\t%.4f%n", pair.getGeneId1(), pair.getGeneId2(), identity, dist);
                    }
                }
            }
        }
        System.out.println("[Phylogeny] Distances saved to: " + outputPrefix + "_ortholog_distances.tsv");
    }

    /**
     * Simple Identity calculation (p-distance)
     */
    private double calculateIdentity(String s1, String s2) {
        int len = Math.min(s1.length(), s2.length());
        int matches = 0;
        for (int i = 0; i < len; i++) {
            if (s1.charAt(i) == s2.charAt(i)) matches++;
        }
        return (double) matches / len;
    }

    /**
     * Jukes-Cantor Distance (d = -3/4 * ln(1 - 4/3 * p))
     */
    private double calculateDistance(String s1, String s2) {
        double identity = calculateIdentity(s1, s2);
        double p = 1.0 - identity;
        if (p >= 0.75) return 3.0; // Saturation limit
        if (p == 0) return 0.0;
        return -0.75 * Math.log(1.0 - (4.0/3.0) * p);
    }
}
