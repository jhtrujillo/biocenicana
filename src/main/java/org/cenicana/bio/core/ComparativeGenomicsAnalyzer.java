package org.cenicana.bio.core;

import org.cenicana.bio.io.CollinearityParser;
import org.cenicana.bio.io.FastaReader;
import org.cenicana.bio.io.GffParser;
import org.cenicana.bio.model.Gene;
import org.cenicana.bio.model.SyntenicBlock;
import org.cenicana.bio.model.SyntenicPair;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Core engine for Comparative Genomics Analysis.
 * Integrates GFF, Collinearity, and FASTA data.
 */
public class ComparativeGenomicsAnalyzer {

    public void runAnalysis(String gff1, String gff2, String collinearity, 
                            String cds1, String cds2, String prot1, String prot2, 
                            String outputTsv) throws IOException {
        
        System.out.println("[CompGen] Loading data...");
        
        GffParser gffParser = new GffParser();
        Map<String, Gene> genes1 = loadGenes(gffParser.parse(gff1));
        Map<String, Gene> genes2 = loadGenes(gffParser.parse(gff2));
        
        CollinearityParser colParser = new CollinearityParser();
        List<SyntenicBlock> blocks = colParser.parse(collinearity);
        
        FastaReader fastaReader = new FastaReader();
        Map<String, String> seqCds1 = cds1 != null ? fastaReader.read(cds1) : Collections.emptyMap();
        Map<String, String> seqCds2 = cds2 != null ? fastaReader.read(cds2) : Collections.emptyMap();
        Map<String, String> seqProt1 = prot1 != null ? fastaReader.read(prot1) : Collections.emptyMap();
        Map<String, String> seqProt2 = prot2 != null ? fastaReader.read(prot2) : Collections.emptyMap();

        System.out.println("[CompGen] Integrating results...");
        
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputTsv))) {
            // Header
            pw.println("Block_ID\tGene1_ID\tChr1\tStart1\tEnd1\tStrand1\tGene2_ID\tChr2\tStart2\tEnd2\tStrand2\tE_Value\tCDS1_Len\tCDS2_Len\tProt1_Len\tProt2_Len");

            for (SyntenicBlock block : blocks) {
                for (SyntenicPair pair : block.getPairs()) {
                    String g1Id = pair.getGeneId1();
                    String g2Id = pair.getGeneId2();
                    
                    Gene g1 = genes1.get(g1Id);
                    Gene g2 = genes2.get(g2Id);
                    
                    if (g1 == null || g2 == null) {
                        // Some IDs might differ between GFF and Collinearity depending on tool settings
                        // Try to find if there's a version suffix or similar (common in CoGe/McScanX)
                        g1 = findFuzzy(genes1, g1Id);
                        g2 = findFuzzy(genes2, g2Id);
                    }

                    String c1 = g1 != null ? g1.getChromosome() : "NA";
                    long s1 = g1 != null ? g1.getStart() : -1;
                    long e1 = g1 != null ? g1.getEnd() : -1;
                    String str1 = g1 != null ? g1.getStrand() : ".";

                    String c2 = g2 != null ? g2.getChromosome() : "NA";
                    long s2 = g2 != null ? g2.getStart() : -1;
                    long e2 = g2 != null ? g2.getEnd() : -1;
                    String str2 = g2 != null ? g2.getStrand() : ".";

                    int cds1Len = seqCds1.getOrDefault(g1Id, "").length();
                    int cds2Len = seqCds2.getOrDefault(g2Id, "").length();
                    int prot1Len = seqProt1.getOrDefault(g1Id, "").length();
                    int prot2Len = seqProt2.getOrDefault(g2Id, "").length();

                    pw.printf(Locale.US, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%e\t%d\t%d\t%d\t%d\n",
                            block.getBlockId(), g1Id, c1, s1, e1, str1, g2Id, c2, s2, e2, str2, 
                            pair.geteValue(), cds1Len, cds2Len, prot1Len, prot2Len);
                }
            }
        }

        System.out.println("[CompGen] Analysis complete. Integrated report saved to: " + outputTsv);
        System.out.println("[CompGen] Summary: " + blocks.size() + " blocks processed.");
    }

    private Map<String, Gene> loadGenes(List<Gene> list) {
        Map<String, Gene> map = new HashMap<>();
        for (Gene g : list) {
            map.put(g.getId(), g);
        }
        return map;
    }

    private Gene findFuzzy(Map<String, Gene> map, String id) {
        if (map.containsKey(id)) return map.get(id);
        // Try common variations like id.1 or transcript:id
        for (String key : map.keySet()) {
            if (key.contains(id) || id.contains(key)) return map.get(key);
        }
        return null;
    }
}
