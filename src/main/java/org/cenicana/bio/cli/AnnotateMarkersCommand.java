package org.cenicana.bio.cli;

import org.cenicana.bio.core.GeneFeature;
import org.cenicana.bio.core.VariantEffectPredictor;
import org.cenicana.bio.io.GeneAnnotationParser;
import org.cenicana.bio.io.AnnotationDashboardGenerator;
import org.cenicana.bio.io.VcfFastReader;
import org.cenicana.bio.io.FastaExtractor;
import org.cenicana.bio.utils.ChromosomeNormalizer;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "annotate", 
         mixinStandardHelpOptions = true,
         version = "1.0",
         description = "Advanced Functional & Population Dashboard for genomic markers.")
public class AnnotateMarkersCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, description = "Input VCF file.")
    private File vcfFile;

    @Option(names = {"-g", "--gff"}, description = "Primary GFF3 file (e.g. 1940).", required = true)
    private File gffFile;

    @Option(names = {"-r", "--ref-genome"}, description = "Reference Genome FASTA file (needed for KASP).")
    private File refGenome;

    @Option(names = {"--gff2"}, description = "Secondary GFF3 file for synteny (e.g. R570).")
    private File gff2File;

    @Option(names = {"-p", "--protein"}, description = "Protein FASTA file.")
    private File proteinFasta;

    @Option(names = {"-c", "--cds"}, description = "CDS FASTA file.")
    private File cdsFasta;

    @Option(names = {"-w", "--window"}, description = "Search window in bp", defaultValue = "5000")
    private long window;

    @Option(names = {"-o", "--output"}, description = "Output HTML file", defaultValue = "genomic_suite.html")
    private String outputFile;

    @Option(names = {"--markers"}, split = ",", description = "Manual marker list (Chr:Pos).")
    private String[] manualMarkers;

    @Override
    public Integer call() throws Exception {
        System.err.println("[BioJava] Loading Primary GFF: " + gffFile.getName());
        Map<String, List<GeneFeature>> index1 = GeneAnnotationParser.parseGenes(gffFile.getAbsolutePath());
        
        Map<String, List<GeneFeature>> index2 = null;
        if (gff2File != null && gff2File.exists()) {
            System.err.println("[BioJava] Loading Secondary GFF (Synteny): " + gff2File.getName());
            index2 = GeneAnnotationParser.parseGenes(gff2File.getAbsolutePath());
        }

        FastaExtractor protExt = (proteinFasta != null && proteinFasta.exists()) ? new FastaExtractor(proteinFasta) : null;
        FastaExtractor cdsExt = (cdsFasta != null && cdsFasta.exists()) ? new FastaExtractor(cdsFasta) : null;
        FastaExtractor genomeExt = (refGenome != null && refGenome.exists()) ? new FastaExtractor(refGenome) : null;
        
        List<AnnotationDashboardGenerator.MarkerMatch> allMatches = new ArrayList<>();
        List<String> samples = new ArrayList<>();
        Map<String, TreeSet<Long>> snpPositions = new HashMap<>();

        if (vcfFile != null && vcfFile.exists()) {
            samples = Arrays.asList(VcfFastReader.getSampleIds(vcfFile.getAbsolutePath()));
            
            // First pass: Index all SNP positions for KASP proximity check
            System.err.println("[BioJava] Indexing SNPs for KASP analysis...");
            for (String[] cols : VcfFastReader.iterateDataBlocks(vcfFile.getAbsolutePath())) {
                String chr = ChromosomeNormalizer.normalize(cols[0]);
                long pos = Long.parseLong(cols[1]);
                snpPositions.computeIfAbsent(chr, k -> new TreeSet<>()).add(pos);
            }

            System.err.println("[BioJava] Processing " + samples.size() + " samples and calculating KASP scores...");
            for (String[] cols : VcfFastReader.iterateDataBlocks(vcfFile.getAbsolutePath())) {
                String chr = ChromosomeNormalizer.normalize(cols[0]);
                long pos = Long.parseLong(cols[1]);
                String markerId = cols[2].equals(".") ? (cols[0] + ":" + cols[1]) : cols[2];
                String ref = cols[3];
                String alt = cols[4];
                
                // Frequency and Genotypes
                int altCount = 0;
                int totalGT = 0;
                String[] gts = new String[samples.size()];
                for (int i = 0; i < samples.size(); i++) {
                    String raw = cols[i+9].split(":")[0];
                    gts[i] = raw;
                    if (raw.contains("1")) altCount++;
                    if (!raw.startsWith(".")) totalGT++;
                }
                double freq = totalGT > 0 ? (double)altCount / totalGT : 0;

                // KASP Calculation
                int kaspScore = 100;
                String flanking = "N/A";
                if (genomeExt != null) {
                    long start = Math.max(1, pos - 50);
                    long end = pos + 50;
                    String seq = genomeExt.extract(chr, start, end);
                    if (seq != null && seq.length() > 5) {
                        int snpPosInSeq = (int)(pos - start);
                        String left = seq.substring(0, snpPosInSeq).toUpperCase();
                        String right = (snpPosInSeq + 1 < seq.length()) ? seq.substring(snpPosInSeq + 1).toUpperCase() : "";
                        flanking = left + "[" + ref + "/" + alt + "]" + right;
                        
                        // Penalty for nearby SNPs (within 30 bp)
                        NavigableSet<Long> nearby = snpPositions.get(chr).subSet(pos - 30, false, pos + 30, false);
                        kaspScore -= (nearby.size() * 20); // -20 per nearby SNP
                        // Penalty for extreme GC (ideal 40-60%)
                        double gc = calculateGC(seq);
                        if (gc < 0.3 || gc > 0.7) kaspScore -= 30;
                        kaspScore = Math.max(0, kaspScore);
                    }
                }

                // Search in Genome 1
                findAndAdd(chr, pos, markerId, ref, alt, freq, gts, kaspScore, flanking, index1, allMatches, protExt, cdsExt, "Primary");
                
                // Search in Genome 2 (Synteny)
                if (index2 != null) {
                    findAndAdd(chr, pos, markerId, ref, alt, freq, gts, kaspScore, flanking, index2, allMatches, null, null, "Secondary");
                }
            }
        }

        System.err.println("[BioJava] Generating Suite: " + outputFile);
        AnnotationDashboardGenerator generator = new AnnotationDashboardGenerator();
        generator.setSamples(samples);
        
        // Collect all unique functional terms from candidate genes
        Set<String> uniqueFunctions = new TreeSet<>();
        for (AnnotationDashboardGenerator.MarkerMatch m : allMatches) {
            uniqueFunctions.addAll(m.gene.getGoTerms());
            uniqueFunctions.addAll(m.gene.getPfamDomains());
            uniqueFunctions.addAll(m.gene.getInterproDomains());
            if (m.gene.getNote() != null && m.gene.getNote().length() > 3) {
                // Add the note as a searchable term if it's short, or its first part
                String note = m.gene.getNote();
                if (note.length() < 50) uniqueFunctions.add(note);
                else uniqueFunctions.add(note.substring(0, 47) + "...");
            }
        }
        generator.setUniqueFunctions(uniqueFunctions);
        
        generator.generate(allMatches, outputFile);
        
        return 0;
    }

    private void findAndAdd(String chr, long pos, String markerId, String ref, String alt, double freq, String[] gts,
                            int kaspScore, String flanking,
                            Map<String, List<GeneFeature>> index, List<AnnotationDashboardGenerator.MarkerMatch> matches,
                            FastaExtractor protExt, FastaExtractor cdsExt, String source) {
        List<GeneFeature> genes = index.get(chr);
        if (genes == null) return;

        long start = pos - window;
        long end = pos + window;

        for (GeneFeature gene : genes) {
            if (gene.overlaps(start, end)) {
                if (protExt != null && gene.getProteinSequence() == null) gene.setProteinSequence(protExt.extract(gene.getId()));
                if (cdsExt != null && gene.getCdsSequence() == null) gene.setCdsSequence(cdsExt.extract(gene.getId()));

                AnnotationDashboardGenerator.MarkerMatch m = new AnnotationDashboardGenerator.MarkerMatch(markerId, chr, pos, gene);
                m.refAllele = ref; m.altAllele = alt; m.alleleFreq = freq; m.genotypes = gts; m.source = source;
                m.kaspScore = kaspScore; m.flankingSeq = flanking;

                if (ref != null && !ref.equals("N") && alt != null && !alt.equals("N")) {
                    VariantEffectPredictor.EffectResult res = VariantEffectPredictor.predict(gene, pos, ref, alt);
                    m.effect = res.type + " (" + res.detail + ")";
                } else {
                    m.effect = "N/A";
                }
                matches.add(m);
            }
            if (gene.getStart() > end) break;
        }
    }

    private double calculateGC(String seq) {
        if (seq == null || seq.isEmpty()) return 0;
        long gc = seq.chars().filter(c -> c == 'G' || c == 'C' || c == 'g' || c == 'c').count();
        return (double) gc / seq.length();
    }
}
