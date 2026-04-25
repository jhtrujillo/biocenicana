package org.cenicana.bio.model;

import java.util.HashMap;
import java.util.Map;

/**
 * Represents a genomic feature (Gene, mRNA, Exon, etc.) from a GFF file.
 */
public class Gene {
    private String id;
    private String chromosome;
    private long start;
    private long end;
    private String strand; // +, -, .
    private String type; // gene, mRNA, exon, etc.
    private Map<String, String> attributes;

    public Gene(String id, String chromosome, long start, long end, String strand, String type) {
        this.id = id;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.type = type;
        this.attributes = new HashMap<>();
    }

    public String getId() { return id; }
    public String getChromosome() { return chromosome; }
    public long getStart() { return start; }
    public long getEnd() { return end; }
    public String getStrand() { return strand; }
    public String getType() { return type; }
    public Map<String, String> getAttributes() { return attributes; }

    /** Returns the first non-empty functional description found in the attributes. */
    public String getDescription() {
        for (String key : new String[]{"description", "Description", "Note", "note", "product", "function", "Dbxref", "Ontology_term"}) {
            String val = attributes.get(key);
            if (val != null && !val.isBlank()) return java.net.URLDecoder.decode(val, java.nio.charset.StandardCharsets.UTF_8);
        }
        return "";
    }

    public void addAttribute(String key, String value) {
        this.attributes.put(key, value);
    }

    @Override
    public String toString() {
        return String.format("Gene[%s, %s:%d-%d (%s)]", id, chromosome, start, end, strand);
    }
}
