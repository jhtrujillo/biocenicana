package org.cenicana.bio.model;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a syntenic block containing multiple collinear gene pairs.
 */
public class SyntenicBlock {
    private String blockId;
    private int numGenes;
    private double score;
    private double eValue;
    private String orientation = "plus";
    private List<SyntenicPair> pairs;
    private boolean hasSV = false;

    public SyntenicBlock(String blockId) {
        this.blockId = blockId;
        this.pairs = new ArrayList<>();
    }

    public String getBlockId() { return blockId; }
    public int getNumGenes() { return numGenes; }
    public void setNumGenes(int numGenes) { this.numGenes = numGenes; }
    public double getScore() { return score; }
    public void setScore(double score) { this.score = score; }
    public double geteValue() { return eValue; }
    public void seteValue(double eValue) { this.eValue = eValue; }
    public String getOrientation() { return orientation; }
    public void setOrientation(String orientation) { this.orientation = orientation; }
    public List<SyntenicPair> getPairs() { return pairs; }
    public boolean hasSV() { return hasSV; }
    public void setHasSV(boolean hasSV) { this.hasSV = hasSV; }

    public void addPair(SyntenicPair pair) {
        this.pairs.add(pair);
    }

    @Override
    public String toString() {
        return String.format("Block[%s, %d genes, score=%.1f]", blockId, numGenes, score);
    }
}
