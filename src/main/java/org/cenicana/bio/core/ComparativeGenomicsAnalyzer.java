package org.cenicana.bio.core;

import org.cenicana.bio.io.AnnotationLoader;
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
                            String outputTsv, String vizOutput,
                            String annotFile1, String annotFile2) throws IOException {
        
        System.out.println("[Phase 1/4] Loading GFF files...");
        GffParser gffParser = new GffParser();
        Map<String, Gene> genes1 = loadGenes(gffParser.parse(gff1));
        Map<String, Gene> genes2 = loadGenes(gffParser.parse(gff2));
        System.out.println("  - Gen1: " + genes1.size() + " genes loaded.");
        System.out.println("  - Gen2: " + genes2.size() + " genes loaded.");

        System.out.println("[Phase 1b/4] Loading Functional Annotations...");
        AnnotationLoader annotLoader = new AnnotationLoader();
        Map<String, String> annot1 = annotLoader.load(annotFile1);
        Map<String, String> annot2 = annotLoader.load(annotFile2);
        System.out.println("  - Annot1: " + annot1.size() + " entries.");
        System.out.println("  - Annot2: " + annot2.size() + " entries.");
        
        System.out.println("[Phase 2/4] Loading Sequence Data (FASTA)...");
        FastaReader fastaReader = new FastaReader();
        Map<String, String> seqCds1 = cds1 != null ? fastaReader.read(cds1) : Collections.emptyMap();
        Map<String, String> seqCds2 = cds2 != null ? fastaReader.read(cds2) : Collections.emptyMap();
        Map<String, String> seqProt1 = prot1 != null ? fastaReader.read(prot1) : Collections.emptyMap();
        Map<String, String> seqProt2 = prot2 != null ? fastaReader.read(prot2) : Collections.emptyMap();
        if (!seqCds1.isEmpty()) System.out.println("  - CDS Gen1: " + seqCds1.size());
        if (!seqCds2.isEmpty()) System.out.println("  - CDS Gen2: " + seqCds2.size());

        System.out.println("[Phase 3/4] Parsing Collinearity and Matching Genes...");
        CollinearityParser colParser = new CollinearityParser();
        List<SyntenicBlock> blocks = colParser.parse(collinearity);
        System.out.println("  - Found " + blocks.size() + " syntenic blocks.");

        Set<String> pairedG1 = new HashSet<>();
        Set<String> pairedG2 = new HashSet<>();

        System.out.println("[Phase 4/4] Integrating and Identifying Orphans...");
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputTsv))) {
            // Header
            pw.println("Block_ID\tStatus\tGene1_ID\tChr1\tStart1\tEnd1\tStrand1\tGene2_ID\tChr2\tStart2\tEnd2\tStrand2\tE_Value\tCDS1_Len\tCDS2_Len\tProt1_Len\tProt2_Len");

            for (SyntenicBlock block : blocks) {
                for (SyntenicPair pair : block.getPairs()) {
                    String g1IdOrig = pair.getGeneId1();
                    String g2IdOrig = pair.getGeneId2();
                    
                    Gene g1 = genes1.get(g1IdOrig);
                    if (g1 == null) g1 = findFuzzy(genes1, g1IdOrig);
                    
                    Gene g2 = genes2.get(g2IdOrig);
                    if (g2 == null) g2 = findFuzzy(genes2, g2IdOrig);
                    
                    if (g1 != null) pairedG1.add(g1.getId());
                    if (g2 != null) pairedG2.add(g2.getId());

                    String g1Id = g1 != null ? g1.getId() : g1IdOrig;
                    String g2Id = g2 != null ? g2.getId() : g2IdOrig;

                    pw.printf(Locale.US, "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%.2e\t%d\t%d\t%d\t%d%n",
                            block.getBlockId(), "Syntenic", g1Id, 
                            g1 != null ? g1.getChromosome() : "NA", g1 != null ? g1.getStart() : -1, g1 != null ? g1.getEnd() : -1, g1 != null ? g1.getStrand() : ".",
                            g2Id,
                            g2 != null ? g2.getChromosome() : "NA", g2 != null ? g2.getStart() : -1, g2 != null ? g2.getEnd() : -1, g2 != null ? g2.getStrand() : ".",
                            pair.geteValue(), 
                            seqCds1.getOrDefault(g1Id, "").length(), seqCds2.getOrDefault(g2Id, "").length(),
                            seqProt1.getOrDefault(g1Id, "").length(), seqProt2.getOrDefault(g2Id, "").length());
                }
            }

            // Identify Orphans in G1
            for (Gene g : genes1.values()) {
                if (!pairedG1.contains(g.getId())) {
                    pw.printf(Locale.US, "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d%n",
                            "None", "Orphan_G1", g.getId(), g.getChromosome(), g.getStart(), g.getEnd(), g.getStrand(),
                            "NA", "NA", -1, -1, ".", "NA",
                            seqCds1.getOrDefault(g.getId(), "").length(), 0, seqProt1.getOrDefault(g.getId(), "").length(), 0);
                }
            }
            // Identify Orphans in G2
            for (Gene g : genes2.values()) {
                if (!pairedG2.contains(g.getId())) {
                    pw.printf(Locale.US, "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d%n",
                            "None", "Orphan_G2", "NA", "NA", -1, -1, ".", 
                            g.getId(), g.getChromosome(), g.getStart(), g.getEnd(), g.getStrand(), "NA",
                            0, seqCds2.getOrDefault(g.getId(), "").length(), 0, seqProt2.getOrDefault(g.getId(), "").length());
                }
            }
        }

        System.out.println("[Done] Integrated report saved to: " + outputTsv);

        if (vizOutput != null) {
            System.out.println("[Viz] Generating interactive visualization...");
            generateVisualization(blocks, genes1, genes2, annot1, annot2, vizOutput);
            System.out.println("[Viz] Visualization saved to: " + vizOutput);
        }
    }

    private void generateVisualization(List<SyntenicBlock> blocks, Map<String, Gene> map1, Map<String, Gene> map2,
                                       Map<String, String> annot1, Map<String, String> annot2,
                                       String outputPath) throws IOException {
        StringBuilder dataJson = new StringBuilder("[");
        boolean first = true;
        for (SyntenicBlock block : blocks) {
            for (SyntenicPair pair : block.getPairs()) {
                Gene g1 = findFuzzy(map1, pair.getGeneId1());
                Gene g2 = findFuzzy(map2, pair.getGeneId2());
                if (g1 != null && g2 != null) {
                    if (!first) dataJson.append(",");
                    // Resolve description: first from dedicated annot map, then from gene attributes
                    String desc1 = annot1.getOrDefault(g1.getId(), g1.getDescription()).replace("\"", "'");
                    String desc2 = annot2.getOrDefault(g2.getId(), g2.getDescription()).replace("\"", "'");
                    dataJson.append(String.format("{\"b\":\"%s\",\"o\":\"%s\",\"g1\":\"%s\",\"c1\":\"%s\",\"s1\":%d,\"e1\":%d,\"f1\":\"%s\",\"g2\":\"%s\",\"c2\":\"%s\",\"s2\":%d,\"e2\":%d,\"f2\":\"%s\"}",
                            block.getBlockId(), block.getOrientation(),
                            g1.getId(), g1.getChromosome(), g1.getStart(), g1.getEnd(), desc1,
                            g2.getId(), g2.getChromosome(), g2.getStart(), g2.getEnd(), desc2));
                    first = false;
                }
            }
        }

        dataJson.append("]");

        // Load HTML template from resources
        String template;
        try (java.io.InputStream is = getClass().getClassLoader().getResourceAsStream("synteny_template.html")) {
            if (is == null) throw new IOException("Template not found: synteny_template.html");
            template = new String(is.readAllBytes(), java.nio.charset.StandardCharsets.UTF_8);
        }
        String html = template.replace("/*DATA_JSON*/", dataJson.toString());

        try (PrintWriter writer = new PrintWriter(new FileWriter(outputPath))) {
            writer.print(html);
        }
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
