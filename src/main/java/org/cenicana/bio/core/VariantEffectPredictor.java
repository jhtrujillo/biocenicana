package org.cenicana.bio.core;

import java.util.*;

/**
 * Predicts the functional effect of a SNP on a protein-coding gene.
 */
public class VariantEffectPredictor {
    
    private static final Map<String, String> GENETIC_CODE = new HashMap<>();
    static {
        String[] codes = {
            "TTT:F", "TTC:F", "TTA:L", "TTG:L", "TCT:S", "TCC:S", "TCA:S", "TCG:S",
            "TAT:Y", "TAC:Y", "TAA:*", "TAG:*", "TGT:C", "TGC:C", "TGA:*", "TGG:W",
            "CTT:L", "CTC:L", "CTA:L", "CTG:L", "CCT:P", "CCC:P", "CCA:P", "CCG:P",
            "CAT:H", "CAC:H", "CAA:Q", "CAG:Q", "CGT:R", "CGC:R", "CGA:R", "CGG:R",
            "ATT:I", "ATC:I", "ATA:I", "ATG:M", "ACT:T", "ACC:T", "ACA:T", "ACG:T",
            "AAT:N", "AAC:N", "AAA:K", "AAG:K", "AGT:S", "AGC:S", "AGA:R", "AGG:R",
            "GTT:V", "GTC:V", "GTA:V", "GTG:V", "GCT:A", "GCC:A", "GCA:A", "GCG:A",
            "GAT:D", "GAC:D", "GAA:E", "GAG:E", "GGT:G", "GGC:G", "GGA:G", "GGG:G"
        };
        for (String c : codes) {
            String[] p = c.split(":");
            GENETIC_CODE.put(p[0], p[1]);
        }
    }

    public static class EffectResult {
        public String type; // Synonymous, Missense, Stop-Gain, Intronic, Non-Coding
        public String detail; // e.g. "V34I"

        public EffectResult(String type, String detail) {
            this.type = type;
            this.detail = detail;
        }
    }

    public static EffectResult predict(GeneFeature gene, long snpPos, String refAllele, String altAllele) {
        if (gene.getCdsSequence() == null || gene.getCdsSequence().isEmpty()) {
            return new EffectResult("Non-Coding", "No CDS available");
        }

        List<GeneFeature.SubFeature> cdsList = new ArrayList<>();
        for (GeneFeature.SubFeature sub : gene.getSubFeatures()) {
            if (sub.type.equalsIgnoreCase("cds")) cdsList.add(sub);
        }

        if (cdsList.isEmpty()) return new EffectResult("Unknown", "No CDS features in GFF");

        // Sort CDS by position
        cdsList.sort(Comparator.comparingLong(s -> s.start));

        // Find which CDS feature contains the SNP
        int cdsOffset = 0;
        GeneFeature.SubFeature targetCds = null;
        for (GeneFeature.SubFeature cds : cdsList) {
            if (snpPos >= cds.start && snpPos <= cds.end) {
                targetCds = cds;
                break;
            }
            cdsOffset += (cds.end - cds.start + 1);
        }

        if (targetCds == null) return new EffectResult("Intronic/UTR", "SNP outside CDS");

        // Calculate relative position in the full CDS string
        int relPos;
        boolean isMinus = "-".equals(gene.getStrand());
        
        if (!isMinus) {
            relPos = cdsOffset + (int)(snpPos - targetCds.start);
        } else {
            // For minus strand, the offset calculation is reversed
            // But usually the CDS sequence we have is already reverse-complemented
            // We need the relative position from the END of the spliced exons
            int totalLen = 0;
            for (GeneFeature.SubFeature cds : cdsList) totalLen += (cds.end - cds.start + 1);
            relPos = totalLen - 1 - (cdsOffset + (int)(snpPos - targetCds.start));
        }

        String cdsSeq = gene.getCdsSequence();
        if (relPos < 0 || relPos >= cdsSeq.length()) return new EffectResult("Error", "Coord out of range");

        // Find codon
        int codonStart = (relPos / 3) * 3;
        if (codonStart + 3 > cdsSeq.length()) return new EffectResult("Error", "Truncated codon");

        String originalCodon = cdsSeq.substring(codonStart, codonStart + 3);
        
        // Mutate codon
        char[] mutatedChars = originalCodon.toCharArray();
        char mutBase = altAllele.charAt(0);
        if (isMinus) mutBase = reverseComplement(mutBase);
        mutatedChars[relPos % 3] = mutBase;
        String mutatedCodon = new String(mutatedChars);

        String aaOrig = GENETIC_CODE.getOrDefault(originalCodon.toUpperCase(), "?");
        String aaMut = GENETIC_CODE.getOrDefault(mutatedCodon.toUpperCase(), "?");

        if (aaOrig.equals(aaMut)) return new EffectResult("Synonymous", aaOrig + (codonStart/3 + 1) + aaMut);
        if (aaMut.equals("*")) return new EffectResult("Stop-Gain", aaOrig + (codonStart/3 + 1) + "*");
        return new EffectResult("Missense", aaOrig + (codonStart/3 + 1) + aaMut);
    }

    private static char reverseComplement(char base) {
        switch (Character.toUpperCase(base)) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: return base;
        }
    }
}
