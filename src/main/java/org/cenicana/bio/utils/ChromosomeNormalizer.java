package org.cenicana.bio.utils;

/**
 * Utility to normalize chromosome names to ensure compatibility between different data sources (VCF, GFF, Fasta).
 */
public class ChromosomeNormalizer {
    
    /**
     * Normalizes a chromosome name:
     * 1. Lowercase
     * 2. Remove 'chr' or 'chromosome' prefixes
     * 3. Remove leading zeros (01 -> 1)
     */
    public static String normalize(String chr) {
        if (chr == null) return "";
        String normalized = chr.toLowerCase().trim();
        
        // Remove common prefixes
        if (normalized.startsWith("chromosome")) {
            normalized = normalized.substring(10);
        } else if (normalized.startsWith("chr")) {
            normalized = normalized.substring(3);
        }
        
        // Remove separators like _ or -
        if (normalized.startsWith("_") || normalized.startsWith("-")) {
            normalized = normalized.substring(1);
        }
        
        // Remove leading zeros (e.g., 01 -> 1)
        while (normalized.length() > 1 && normalized.startsWith("0")) {
            normalized = normalized.substring(1);
        }
        
        return normalized;
    }
}
