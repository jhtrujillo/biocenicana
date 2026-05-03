package org.cenicana.bio.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Represents a genomic feature (usually a gene) extracted from a GFF file.
 * Stores spatial coordinates and functional metadata (GO, PFAM, etc.)
 */
public class GeneFeature {
    private String chromosome;
    private long start;
    private long end;
    private String strand;
    private String id;
    private String name;
    private String note;
    
    // Sub-features like Exons and CDS for mapping
    private List<SubFeature> subFeatures = new ArrayList<>();
    
    // Functional terms extracted from Column 9
    private Set<String> goTerms = new LinkedHashSet<>();
    private Set<String> pfamDomains = new LinkedHashSet<>();
    private Set<String> interproDomains = new LinkedHashSet<>();
    private Map<String, String> otherAttributes = new HashMap<>();
    
    private String proteinSequence;
    private String cdsSequence;

    public static class SubFeature {
        public String type;
        public long start;
        public long end;
        public SubFeature(String type, long start, long end) {
            this.type = type;
            this.start = Math.min(start, end);
            this.end = Math.max(start, end);
        }
    }

    public GeneFeature(String chromosome, long start, long end, String strand) {
        this.chromosome = chromosome;
        this.start = Math.min(start, end);
        this.end = Math.max(start, end);
        this.strand = strand;
    }

    // --- Getters and Setters ---
    public String getChromosome() { return chromosome; }
    public long getStart() { return start; }
    public long getEnd() { return end; }
    public String getStrand() { return strand; }
    public String getId() { return id; }
    public void setId(String id) { this.id = id; }
    public String getName() { return name; }
    public void setName(String name) { this.name = name; }
    public String getNote() { return note; }
    public void setNote(String note) { this.note = note; }
    public List<SubFeature> getSubFeatures() { return subFeatures; }
    public void addSubFeature(String type, long start, long end) {
        subFeatures.add(new SubFeature(type, start, end));
    }
    public Set<String> getGoTerms() { return goTerms; }
    public Set<String> getPfamDomains() { return pfamDomains; }
    public Set<String> getInterproDomains() { return interproDomains; }
    public Map<String, String> getOtherAttributes() { return otherAttributes; }

    public String getProteinSequence() { return proteinSequence; }
    public void setProteinSequence(String proteinSequence) { this.proteinSequence = proteinSequence; }
    public String getCdsSequence() { return cdsSequence; }
    public void setCdsSequence(String cdsSequence) { this.cdsSequence = cdsSequence; }

    public void addGoTerm(String term) { if (term != null && !term.isEmpty()) goTerms.add(term); }
    public void addPfamDomain(String domain) { if (domain != null && !domain.isEmpty()) pfamDomains.add(domain); }
    public void addInterproDomain(String domain) { if (domain != null && !domain.isEmpty()) interproDomains.add(domain); }

    /**
     * Checks if this gene overlaps with a given genomic range.
     */
    public boolean overlaps(long rangeStart, long rangeEnd) {
        return this.start <= rangeEnd && this.end >= rangeStart;
    }

    /**
     * Calculates the minimum distance from a point to this gene.
     * Returns 0 if the point is inside the gene.
     */
    public long distanceTo(long position) {
        if (position >= start && position <= end) return 0;
        return Math.min(Math.abs(position - start), Math.abs(position - end));
    }
}
