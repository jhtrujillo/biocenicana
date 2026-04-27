package org.cenicana.bio.core;

import org.cenicana.bio.io.AnnotationLoader;
import org.cenicana.bio.io.CollinearityParser;
import org.cenicana.bio.io.FastaReader;
import org.cenicana.bio.io.GffParser;
import org.cenicana.bio.io.SvParser;
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
                            String kaksFile, String exportOrthologs, String svFile) throws IOException {
        
        System.out.println("[Phase 1/4] Loading GFF files...");
        GffParser gffParser = new GffParser();
        Map<String, Gene> genes1 = loadGenes(gffParser.parse(gff1));
        Map<String, Gene> genes2 = loadGenes(gffParser.parse(gff2));
        System.out.println("  - Gen1: " + genes1.size() + " genes loaded.");
        System.out.println("  - Gen2: " + genes2.size() + " genes loaded.");

        System.out.println("[Phase 1b/4] Loading Functional Annotations...");
        AnnotationLoader annotLoader = new AnnotationLoader();
        String aFile1 = (annotFile1 != null && !annotFile1.isBlank()) ? annotFile1 : gff1;
        String aFile2 = (annotFile2 != null && !annotFile2.isBlank()) ? annotFile2 : gff2;

        Map<String, String> annot1 = annotLoader.load(aFile1);
        Map<String, String> annot2 = annotLoader.load(aFile2);
        System.out.println("  - Annot1: " + annot1.size() + " entries.");
        System.out.println("  - Annot2: " + annot2.size() + " entries.");

        System.out.println("[Phase 1d/4] Loading GO Terms...");
        Map<String, List<String>> go1 = annotLoader.loadGoTerms(aFile1);
        Map<String, List<String>> go2 = annotLoader.loadGoTerms(aFile2);
        System.out.println("  - GO1: " + go1.size() + " genes with GO.");
        System.out.println("  - GO2: " + go2.size() + " genes with GO.");

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

        Map<String, Gene> base1 = createBaseIdMap(genes1);
        Map<String, Gene> base2 = createBaseIdMap(genes2);

        // --- KA/KS CALCULATION / LOADING ---
        Map<String, double[]> kaksData;
        if (kaksFile == null || kaksFile.isBlank()) {
            if (!seqCds1.isEmpty() && !seqCds2.isEmpty()) {
                System.out.println("[Phase 3b/4] Calculating Ka/Ks Selection Pressure (Parallel)...");
                kaksData = calculateKaksParallel(blocks, genes1, genes2, base1, base2, seqCds1, seqCds2);
                saveKaksReport(kaksData, "results/auto_kaks_report.tsv");
            } else {
                kaksData = Collections.emptyMap();
            }
        } else {
            System.out.println("[Phase 1c/4] Loading existing Ka/Ks Data...");
            kaksData = loadKaksFile(kaksFile);
        }

        System.out.println("[Phase 4/4] Integrating and Identifying Orphans...");
        Set<String> pairedG1 = new HashSet<>();
        Set<String> pairedG2 = new HashSet<>();
        
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputTsv))) {
            // Header
            pw.println("Block_ID\tStatus\tGene1_ID\tChr1\tStart1\tEnd1\tStrand1\tGene2_ID\tChr2\tStart2\tEnd2\tStrand2\tE_Value\tCDS1_Len\tCDS2_Len\tProt1_Len\tProt2_Len");

            for (SyntenicBlock block : blocks) {
                for (SyntenicPair pair : block.getPairs()) {
                    String g1IdOrig = pair.getGeneId1();
                    String g2IdOrig = pair.getGeneId2();
                    
                    Gene g1 = fastFind(g1IdOrig, genes1, base1);
                    Gene g2 = fastFind(g2IdOrig, genes2, base2);
                    
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

        // --- ORTHOLOG EXPORT & GO ENRICHMENT ---
        if (exportOrthologs != null && !exportOrthologs.isBlank()) {
            exportOrthologsForPhylogeny(blocks, genes1, genes2, base1, base2, seqCds1, seqCds2, go1, go2, exportOrthologs);
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
                        Gene g1 = fastFind(pair.getGeneId1(), genes1, base1);
                        Gene g2 = fastFind(pair.getGeneId2(), genes2, base2);
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

            // Phase 4d: Structural Variant (SV) Intersection
            if (svFile != null && !svFile.isBlank()) {
                System.out.println("[Phase 4d/5] Intersecting syntenic blocks with Structural Variants (SVs)...");
                SvParser svParser = new SvParser();
                List<SvParser.SvRegion> svs = svParser.parse(svFile);
                System.out.println("  - Loaded " + svs.size() + " SV regions.");
                
                for (SyntenicBlock block : blocks) {
                    // Get bounding box for the block on both genomes
                    String chr1 = null, chr2 = null;
                    long min1 = Long.MAX_VALUE, max1 = 0;
                    long min2 = Long.MAX_VALUE, max2 = 0;
                    
                    for (SyntenicPair pair : block.getPairs()) {
                        Gene g1 = fastFind(pair.getGeneId1(), genes1, base1);
                        Gene g2 = fastFind(pair.getGeneId2(), genes2, base2);
                        if (g1 != null) {
                            chr1 = g1.getChromosome();
                            min1 = Math.min(min1, g1.getStart());
                            max1 = Math.max(max1, g1.getEnd());
                        }
                        if (g2 != null) {
                            chr2 = g2.getChromosome();
                            min2 = Math.min(min2, g2.getStart());
                            max2 = Math.max(max2, g2.getEnd());
                        }
                    }
                    
                    for (SvParser.SvRegion sv : svs) {
                        boolean match1 = chr1 != null && normalizeChromosome(chr1).equals(normalizeChromosome(sv.chromosome)) 
                                         && sv.start <= max1 && sv.end >= min1;
                        boolean match2 = chr2 != null && normalizeChromosome(chr2).equals(normalizeChromosome(sv.chromosome)) 
                                         && sv.start <= max2 && sv.end >= min2;
                        
                        if (match1 || match2) {
                            block.setHasSV(true);
                            break;
                        }
                    }
                }
            }

            System.out.println("[Viz] Generating interactive visualization...");
            String n1 = gff1.toLowerCase().contains("1940") ? "CC 01-1940" : (gff1.toLowerCase().contains("r570") ? "R570" : "Genome 1");
            String n2 = gff2.toLowerCase().contains("1940") ? "CC 01-1940" : (gff2.toLowerCase().contains("r570") ? "R570" : "Genome 2");
            generateVisualization(blocks, genes1, genes2, base1, base2, annot1, annot2, go1, go2, blockDiv, kaksData, n1, n2, vizOutput);
            System.out.println("[Viz] Visualization saved to: " + vizOutput);
        }

        // Ortholog export already handled in Phase 4
    }

    // Deprecated old method, unified logic in exportOrthologsForPhylogeny

    private void generateVisualization(List<SyntenicBlock> blocks, Map<String, Gene> map1, Map<String, Gene> map2,
                                       Map<String, Gene> base1, Map<String, Gene> base2,
                                       Map<String, String> annot1, Map<String, String> annot2,
                                       Map<String, List<String>> go1, Map<String, List<String>> go2,
                                       Map<String, Double> blockDiv,
                                       Map<String, double[]> kaksData,
                                       String n1, String n2,
                                       String outputPath) throws IOException {
        GoEnrichmentCalculator goCalc = new GoEnrichmentCalculator();
        
        // --- PARALLEL GO ENRICHMENT FOR ALL BLOCKS ---
        System.out.println("  - Calculating functional enrichment for " + blocks.size() + " blocks in parallel...");
        Map<String, List<GoEnrichmentCalculator.EnrichmentResult>> blockEnrichment = new java.util.concurrent.ConcurrentHashMap<>();
        
        blocks.parallelStream().forEach(block -> {
            Map<String, List<String>> studyGo = new HashMap<>();
            for (SyntenicPair pair : block.getPairs()) {
                Gene g1 = fastFind(pair.getGeneId1(), map1, base1);
                if (g1 != null && go1.containsKey(g1.getId())) studyGo.put(g1.getId(), go1.get(g1.getId()));
                Gene g2 = fastFind(pair.getGeneId2(), map2, base2);
                if (g2 != null && go2.containsKey(g2.getId())) studyGo.put(g2.getId(), go2.get(g2.getId()));
            }
            if (!studyGo.isEmpty()) {
                Map<String, List<String>> bgGo = new HashMap<>(go1);
                bgGo.putAll(go2);
                List<GoEnrichmentCalculator.EnrichmentResult> enriched = goCalc.calculate(studyGo, bgGo);
                blockEnrichment.put(block.getBlockId(), enriched);
            }
        });

        StringBuilder dataJson = new StringBuilder("[");
        boolean first = true;
        for (SyntenicBlock block : blocks) {
            List<GoEnrichmentCalculator.EnrichmentResult> enriched = blockEnrichment.getOrDefault(block.getBlockId(), Collections.emptyList());
            StringBuilder goJson = new StringBuilder("[");
            for (int i = 0; i < Math.min(enriched.size(), 5); i++) {
                GoEnrichmentCalculator.EnrichmentResult r = enriched.get(i);
                if (i > 0) goJson.append(",");
                goJson.append(String.format(java.util.Locale.US, "{\"id\":\"%s\",\"p\":%.4f,\"c\":%d}", r.goId, r.pValue, r.studyCount));
            }
            goJson.append("]");

            for (SyntenicPair pair : block.getPairs()) {
                Gene g1 = fastFind(pair.getGeneId1(), map1, base1);
                Gene g2 = fastFind(pair.getGeneId2(), map2, base2);
                if (g1 != null && g2 != null) {
                    if (!first) dataJson.append(",");
                    // Resolve description: first from dedicated annot map, then from gene attributes
                    String desc1 = annot1.getOrDefault(g1.getId(), g1.getDescription()).replace("\"", "'");
                    String desc2 = annot2.getOrDefault(g2.getId(), g2.getDescription()).replace("\"", "'");
                    double div = blockDiv.getOrDefault(block.getBlockId(), 0.0);
                    double[] kk = kaksData.getOrDefault(g1.getId() + ":" + g2.getId(), null);
                    String kkJson = kk == null ? "\"kk\":null" : String.format(java.util.Locale.US, "\"kk\":{\"ka\":%.4f,\"ks\":%.4f,\"r\":%.4f}", kk[0], kk[1], kk[2]);
                    
                    dataJson.append(String.format(java.util.Locale.US, "{\"b\":\"%s\",\"o\":\"%s\",\"div\":%.2f,%s,\"sv\":%b,\"go\":%s,\"g1\":\"%s\",\"c1\":\"%s\",\"s1\":%d,\"e1\":%d,\"f1\":\"%s\",\"g2\":\"%s\",\"c2\":\"%s\",\"s2\":%d,\"e2\":%d,\"f2\":\"%s\"}",
                            block.getBlockId(), block.getOrientation(), div, kkJson, block.hasSV(), goJson.toString(),
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

        String html = template.replace("/*DATA_JSON*/", dataJson.toString())
                              .replace("/*G1_NAME*/", n1)
                              .replace("/*G2_NAME*/", n2);

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

    private Map<String, double[]> calculateKaksParallel(List<SyntenicBlock> blocks, Map<String, Gene> map1, Map<String, Gene> map2,
                                                       Map<String, Gene> base1, Map<String, Gene> base2,
                                                       Map<String, String> seqs1, Map<String, String> seqs2) {
        System.out.println("  - Aligning and calculating Ka/Ks for all syntenic pairs...");
        KaKsCalculator calculator = new KaKsCalculator();
        Map<String, double[]> results = new java.util.concurrent.ConcurrentHashMap<>();

        blocks.parallelStream().forEach(block -> {
            for (SyntenicPair pair : block.getPairs()) {
                Gene g1 = fastFind(pair.getGeneId1(), map1, base1);
                Gene g2 = fastFind(pair.getGeneId2(), map2, base2);
                if (g1 != null && g2 != null) {
                    String s1 = seqs1.get(g1.getId());
                    String s2 = seqs2.get(g2.getId());
                    if (s1 != null && s2 != null && !s1.isBlank() && !s2.isBlank()) {
                        try {
                            String[] aligned = calculator.align(s1, s2);
                            KaKsCalculator.Result res = calculator.calculate(aligned[0], aligned[1]);
                            results.put(pair.getGeneId1() + ":" + pair.getGeneId2(), new double[]{res.ka, res.ks, res.ratio});
                        } catch (Exception ignored) {}
                    }
                }
            }
        });
        System.out.println("  - Completed Ka/Ks for " + results.size() + " pairs.");
        return results;
    }

    private void saveKaksReport(Map<String, double[]> data, String path) throws IOException {
        java.nio.file.Files.createDirectories(java.nio.file.Paths.get("results"));
        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            pw.println("#Gene1\tGene2\tKa\tKs\tKa/Ks");
            for (Map.Entry<String, double[]> entry : data.entrySet()) {
                String[] ids = entry.getKey().split(":");
                double[] vals = entry.getValue();
                pw.printf(Locale.US, "%s\t%s\t%.6f\t%.6f\t%.6f%n", ids[0], ids[1], vals[0], vals[1], vals[2]);
            }
        }
    }

    private Map<String, double[]> loadKaksFile(String path) {
        Map<String, double[]> data = new HashMap<>();
        try (java.io.BufferedReader br = java.nio.file.Files.newBufferedReader(java.nio.file.Paths.get(path))) {
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
                        data.put(key, new double[]{ka, ks, ratio});
                    } catch (Exception ignored) {}
                }
            }
        } catch (IOException e) {
            System.err.println("  - Error loading Ka/Ks file: " + e.getMessage());
        }
        return data;
    }

    private Map<String, Gene> createBaseIdMap(Map<String, Gene> originalMap) {
        Map<String, Gene> baseMap = new HashMap<>(originalMap.size());
        for (Gene g : originalMap.values()) {
            String stripped = stripPrefix(g.getId());
            baseMap.putIfAbsent(stripped, g);
            // Also index by ID without version suffix if possible (e.g., Soffi.1 -> Soffi)
            int lastDot = stripped.lastIndexOf('.');
            if (lastDot > 0) {
                baseMap.putIfAbsent(stripped.substring(0, lastDot), g);
            }
        }
        return baseMap;
    }

    private String stripPrefix(String id) {
        if (id == null) return "";
        return id.replaceFirst("(?i)^(gene:|mrna:|transcript:)", "");
    }

    private Gene fastFind(String id, Map<String, Gene> fullMap, Map<String, Gene> baseMap) {
        if (fullMap.containsKey(id)) return fullMap.get(id);
        String stripped = stripPrefix(id);
        if (fullMap.containsKey(stripped)) return fullMap.get(stripped);
        if (baseMap.containsKey(stripped)) return baseMap.get(stripped);
        
        // Try removing suffix
        int lastDot = stripped.lastIndexOf('.');
        if (lastDot > 0) {
            String noVer = stripped.substring(0, lastDot);
            if (baseMap.containsKey(noVer)) return baseMap.get(noVer);
        }
        return null;
    }
    
    private String normalizeChromosome(String chr) {
        if (chr == null) return "";
        // Remove common prefixes ignoring case
        String s = chr.replaceAll("(?i)^(chromosome|chr|contig|scaffold)_?", "");
        // Keep only alphanumeric characters to handle variants like "1A", "01", etc.
        s = s.replaceAll("[^a-zA-Z0-9]", "");
        // Remove leading zeros
        s = s.replaceFirst("^0+(?!$)", "");
        return s.toLowerCase();
    }

    private void exportOrthologsForPhylogeny(List<SyntenicBlock> blocks, Map<String, Gene> map1, Map<String, Gene> map2,
                                            Map<String, Gene> base1, Map<String, Gene> base2,
                                            Map<String, String> seqs1, Map<String, String> seqs2,
                                            Map<String, List<String>> go1, Map<String, List<String>> go2,
                                            String exportDir) throws IOException {
        System.out.println("[Phylogeny] Exporting 1:1 orthologs & Super-Matrix to: " + exportDir);
        java.nio.file.Files.createDirectories(java.nio.file.Paths.get(exportDir));
        
        KaKsCalculator aligner = new KaKsCalculator();
        List<String[]> pairsToAlign = new ArrayList<>();
        Map<String, List<String>> orthologGoSet = new HashMap<>();

        for (SyntenicBlock block : blocks) {
            for (SyntenicPair pair : block.getPairs()) {
                Gene g1 = fastFind(pair.getGeneId1(), map1, base1);
                Gene g2 = fastFind(pair.getGeneId2(), map2, base2);
                if (g1 != null && g2 != null) {
                    String s1 = seqs1.get(g1.getId());
                    String s2 = seqs2.get(g2.getId());
                    if (s1 != null && s2 != null && s1.length() > 300) {
                        String outName = String.format("%s_%s.fasta", g1.getId(), g2.getId());
                        try (PrintWriter pw = new PrintWriter(new FileWriter(new java.io.File(exportDir, outName)))) {
                            pw.printf(">%s%n%s%n", g1.getId(), s1);
                            pw.printf(">%s%n%s%n", g2.getId(), s2);
                        }
                        pairsToAlign.add(new String[]{s1, s2});
                        if (go1.containsKey(g1.getId())) orthologGoSet.put(g1.getId(), go1.get(g1.getId()));
                        if (go2.containsKey(g2.getId())) orthologGoSet.put(g2.getId(), go2.get(g2.getId()));
                    }
                }
            }
            if (pairsToAlign.size() > 5000) break;
        }
        System.out.println("  - Exported " + pairsToAlign.size() + " individual ortholog fasta files.");

        // --- SUPER-MATRIX ---
        if (!pairsToAlign.isEmpty()) {
            System.out.println("  - Aligning orthologs for Super-Matrix...");
            List<String[]> alignedPairs = pairsToAlign.parallelStream().map(p -> aligner.align(p[0], p[1])).collect(java.util.stream.Collectors.toList());
            StringBuilder super1 = new StringBuilder(), super2 = new StringBuilder();
            for (String[] al : alignedPairs) { super1.append(al[0]); super2.append(al[1]); }
            String smPath = new java.io.File(exportDir, "supermatrix_orthologs.fasta").getAbsolutePath();
            try (PrintWriter pw = new PrintWriter(new FileWriter(smPath))) {
                pw.println(">CC_01_1940_SuperMatrix\n" + super1.toString());
                pw.println(">R570_SuperMatrix\n" + super2.toString());
            }
            System.out.println("  - Super-Matrix saved to: " + smPath);
        }

        // --- GO ENRICHMENT ---
        if (!orthologGoSet.isEmpty()) {
            System.out.println("[Functional] Calculating GO Enrichment for orthologs...");
            GoEnrichmentCalculator goCalc = new GoEnrichmentCalculator();
            Map<String, List<String>> bg = new HashMap<>(go1); bg.putAll(go2);
            List<GoEnrichmentCalculator.EnrichmentResult> res = goCalc.calculate(orthologGoSet, bg);
            String goOut = new java.io.File(exportDir, "ortholog_go_enrichment.tsv").getAbsolutePath();
            try (PrintWriter pw = new PrintWriter(new FileWriter(goOut))) {
                pw.println("GO_Term\tCategory\tDescription\tStudy_Count\tBackground_Count\tpValue\tAdjusted_pValue");
                for (GoEnrichmentCalculator.EnrichmentResult r : res) {
                    if (r.getPValue() < 0.05) {
                        pw.printf(Locale.US, "%s\t%s\t%s\t%d\t%d\t%.2e\t%.2e%n", r.getGoId(), r.getCategory(), r.getDescription(), r.getStudyCount(), r.getBgCount(), r.getPValue(), r.getAdjPValue());
                    }
                }
            }
            System.out.println("  - GO Enrichment report saved to: " + goOut);
        }
    }

    private static class VcfRegion {
        String blockId;
        long start, end, count;
        VcfRegion(String id, long s, long e) {
            this.blockId = id; this.start = s; this.end = e; this.count = 0;
        }
    }
}
