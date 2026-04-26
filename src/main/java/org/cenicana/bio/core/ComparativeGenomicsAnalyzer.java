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
                            String annotFile1, String annotFile2, String vcfFile,
                            String kaksFile) throws IOException {
        
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

        System.out.println("[Phase 1c/4] Loading Ka/Ks Selection Data...");
        Map<String, double[]> kaksData = new HashMap<>(); // Key: G1:G2, Value: [Ka, Ks, Ka/Ks]
        if (kaksFile != null && !kaksFile.isBlank()) {
            try (java.io.BufferedReader br = java.nio.file.Files.newBufferedReader(java.nio.file.Paths.get(kaksFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.trim().isEmpty() || line.startsWith("#")) continue;
                    String[] p = line.split("\t");
                    if (p.length >= 5) {
                        try {
                            String key = p[0] + ":" + p[1];
                            double ka = Double.parseDouble(p[2]);
                            double ks = Double.parseDouble(p[3]);
                            double ratio = Double.parseDouble(p[4]);
                            kaksData.put(key, new double[]{ka, ks, ratio});
                        } catch (Exception ignored) {}
                    }
                }
            }
            System.out.println("  - Loaded Ka/Ks for " + kaksData.size() + " gene pairs.");
        }
        
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

        Map<String, Double> blockDiv = new HashMap<>();
        if (vizOutput != null) {
            if (vcfFile != null && !vcfFile.isBlank()) {
                System.out.println("[Phase 4b/5] Calculando abundancia de marcadores en regiones colineales (VCF)...");
                // --- HIGH PERFORMANCE VCF PARSING ---
                System.out.println("[Phase 4b/5] Indexing blocks by chromosome for fast VCF lookup...");
                // Index: Map<NormalizedChr, List<RegionRecord>>
                Map<String, List<VcfRegion>> chrIndex = new HashMap<>();
                for (SyntenicBlock block : blocks) {
                    long minStart1 = Long.MAX_VALUE, maxEnd1 = 0;
                    long minStart2 = Long.MAX_VALUE, maxEnd2 = 0;
                    String chr1 = null, chr2 = null;
                    for (SyntenicPair pair : block.getPairs()) {
                        Gene g1 = findFuzzy(genes1, pair.getGeneId1());
                        Gene g2 = findFuzzy(genes2, pair.getGeneId2());
                        if (g1 != null) {
                            chr1 = g1.getChromosome();
                            minStart1 = Math.min(minStart1, g1.getStart());
                            maxEnd1 = Math.max(maxEnd1, g1.getEnd());
                        }
                        if (g2 != null) {
                            chr2 = g2.getChromosome();
                            minStart2 = Math.min(minStart2, g2.getStart());
                            maxEnd2 = Math.max(maxEnd2, g2.getEnd());
                        }
                    }
                    if (chr1 != null) {
                        String nChr = normalizeChromosome(chr1);
                        chrIndex.computeIfAbsent(nChr, k -> new ArrayList<>()).add(new VcfRegion(block.getBlockId(), minStart1, maxEnd1));
                    }
                    if (chr2 != null) {
                        String nChr = normalizeChromosome(chr2);
                        chrIndex.computeIfAbsent(nChr, k -> new ArrayList<>()).add(new VcfRegion(block.getBlockId(), minStart2, maxEnd2));
                    }
                }

                System.out.println("  - Starting high-speed VCF scan...");
                long snpTotal = 0;
                try (java.io.BufferedReader br = java.nio.file.Files.newBufferedReader(java.nio.file.Paths.get(vcfFile))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        if (line.isEmpty() || line.charAt(0) == '#') continue;
                        
                        // Ultra-fast parsing of first two columns without full split
                        int tab1 = line.indexOf('\t');
                        if (tab1 < 1) continue;
                        int tab2 = line.indexOf('\t', tab1 + 1);
                        if (tab2 < 1) continue;
                        
                        String vcfChrRaw = line.substring(0, tab1);
                        String vcfChr = normalizeChromosome(vcfChrRaw);
                        
                        List<VcfRegion> regions = chrIndex.get(vcfChr);
                        if (regions != null) {
                            long pos = Long.parseLong(line.substring(tab1 + 1, tab2));
                            for (VcfRegion reg : regions) {
                                if (pos >= reg.start && pos <= reg.end) {
                                    reg.count++;
                                }
                            }
                        }
                        snpTotal++;
                        if (snpTotal % 500000 == 0) System.out.print(".");
                    }
                }
                System.out.println("\n  - VCF scan complete. Calculating final density...");

                // Aggregating counts back to blockDiv
                for (List<VcfRegion> regions : chrIndex.values()) {
                    for (VcfRegion reg : regions) {
                        if (reg.count > 0) {
                            double sizeKb = Math.max(1, (reg.end - reg.start) / 1000.0);
                            double density = reg.count / sizeKb;
                            blockDiv.put(reg.blockId, blockDiv.getOrDefault(reg.blockId, 0.0) + density);
                        }
                    }
                }
            }

            System.out.println("[Viz] Generating interactive visualization...");
            generateVisualization(blocks, genes1, genes2, annot1, annot2, blockDiv, kaksData, vizOutput);
            System.out.println("[Viz] Visualization saved to: " + vizOutput);
        }
    }

    private void generateVisualization(List<SyntenicBlock> blocks, Map<String, Gene> map1, Map<String, Gene> map2,
                                       Map<String, String> annot1, Map<String, String> annot2,
                                       Map<String, Double> blockDiv,
                                       Map<String, double[]> kaksData,
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
                    double div = blockDiv.getOrDefault(block.getBlockId(), 0.0);
                    double[] kk = kaksData.getOrDefault(g1.getId() + ":" + g2.getId(), null);
                    String kkJson = kk == null ? "\"kk\":null" : String.format(java.util.Locale.US, "\"kk\":{\"ka\":%.4f,\"ks\":%.4f,\"r\":%.4f}", kk[0], kk[1], kk[2]);
                    
                    dataJson.append(String.format(java.util.Locale.US, "{\"b\":\"%s\",\"o\":\"%s\",\"div\":%.2f,%s,\"g1\":\"%s\",\"c1\":\"%s\",\"s1\":%d,\"e1\":%d,\"f1\":\"%s\",\"g2\":\"%s\",\"c2\":\"%s\",\"s2\":%d,\"e2\":%d,\"f2\":\"%s\"}",
                            block.getBlockId(), block.getOrientation(), div, kkJson,
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
    
    private String normalizeChromosome(String chr) {
        if (chr == null) return "";
        String s = chr.toLowerCase().replaceAll("[^a-z0-9]", "");
        s = s.replace("1940", "").replace("r570", "").replace("cc01", "");
        s = s.replace("chromosome", "").replace("chr", "").replace("contig", "").replace("scaffold", "");
        s = s.replaceFirst("^0+(?!$)", "");
        return s;
    }

    private static class VcfRegion {
        String blockId;
        long start, end, count;
        VcfRegion(String id, long s, long e) {
            this.blockId = id; this.start = s; this.end = e; this.count = 0;
        }
    }
}
