package org.cenicana.bio.core;

import java.util.*;

/**
 * Implements the Neighbor-Joining (NJ) algorithm to build phylogenetic trees
 * from a genetic distance matrix.
 */
public class PhylogenyTreeBuilder {

    public String buildNewick(float[][] matrix, String[] labels) {
        int n = labels.length;
        List<Node> nodes = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            nodes.add(new Node(labels[i]));
        }

        double[][] d = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                d[i][j] = matrix[i][j];
            }
        }

        while (nodes.size() > 2) {
            int size = nodes.size();
            double[] r = new double[size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    r[i] += d[i][j];
                }
            }

            // Find pair (i, j) that minimizes Q
            int bestI = 0, bestJ = 1;
            double minQ = Double.MAX_VALUE;
            for (int i = 0; i < size; i++) {
                for (int j = i + 1; j < size; j++) {
                    double q = (size - 2) * d[i][j] - r[i] - r[j];
                    if (q < minQ) {
                        minQ = q;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }

            // Calculate branch lengths
            double dIJ = d[bestI][bestJ];
            double li = 0.5 * dIJ + (r[bestI] - r[bestJ]) / (2 * (size - 2));
            double lj = dIJ - li;

            // Create new node
            Node newNode = new Node("(" + nodes.get(bestI).label + ":" + String.format(Locale.US, "%.6f", Math.max(0, li)) + "," + 
                                          nodes.get(bestJ).label + ":" + String.format(Locale.US, "%.6f", Math.max(0, lj)) + ")");
            
            // New distances to all other nodes
            double[][] newD = new double[size - 1][size - 1];
            List<Node> nextNodes = new ArrayList<>();
            int nextIdx = 0;
            
            for (int m = 0; m < size; m++) {
                if (m != bestI && m != bestJ) {
                    nextNodes.add(nodes.get(m));
                    nextIdx++;
                }
            }
            nextNodes.add(newNode);

            for (int i = 0; i < size - 1; i++) {
                for (int j = 0; j < size - 1; j++) {
                    if (i < size - 2 && j < size - 2) {
                        // Old distances
                        int oldI = findOldIdx(nodes, nextNodes.get(i), bestI, bestJ);
                        int oldJ = findOldIdx(nodes, nextNodes.get(j), bestI, bestJ);
                        newD[i][j] = d[oldI][oldJ];
                    } else if (i == size - 2 && j < size - 2) {
                        // Distance to new node
                        int m = findOldIdx(nodes, nextNodes.get(j), bestI, bestJ);
                        newD[i][j] = 0.5 * (d[bestI][m] + d[bestJ][m] - dIJ);
                        newD[j][i] = newD[i][j];
                    }
                }
            }

            d = newD;
            nodes = nextNodes;
        }

        // Final two nodes
        return "(" + nodes.get(0).label + ":" + String.format(Locale.US, "%.6f", Math.max(0, d[0][1])) + "," + 
                     nodes.get(1).label + ":0.000000);" ;
    }

    private int findOldIdx(List<Node> oldNodes, Node target, int skip1, int skip2) {
        for (int i = 0; i < oldNodes.size(); i++) {
            if (oldNodes.get(i) == target) return i;
        }
        return -1;
    }

    private static class Node {
        String label;
        Node(String label) { this.label = label; }
    }
}
