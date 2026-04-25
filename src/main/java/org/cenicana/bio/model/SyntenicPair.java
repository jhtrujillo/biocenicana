package org.cenicana.bio.model;

/**
 * Represents a pair of syntenic/collinear genes.
 */
public class SyntenicPair {
    private String geneId1;
    private String geneId2;
    private double eValue;

    public SyntenicPair(String geneId1, String geneId2, double eValue) {
        this.geneId1 = geneId1;
        this.geneId2 = geneId2;
        this.eValue = eValue;
    }

    public String getGeneId1() { return geneId1; }
    public String getGeneId2() { return geneId2; }
    public double geteValue() { return eValue; }

    @Override
    public String toString() {
        return String.format("%s <=> %s (E=%e)", geneId1, geneId2, eValue);
    }
}
