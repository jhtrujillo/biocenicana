package org.cenicana.bio.core;

import org.cenicana.bio.io.VcfFastReader;
import org.ejml.simple.SimpleMatrix;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import org.ejml.data.DMatrixRMaj;

import java.io.*;
import java.util.*;

/**
 * Population Structure Analyzer using Principal Component Analysis (PCA).
 */
public class PopulationStructureAnalyzer {

    public static class PcaResult {
        public String[] sampleNames;
        public double[][] pcMatrix; // [numSamples][numPCs]
        public double[] explainedVariance;
        public double[] eigenvalues;
        public int[] clusterAssignments; // Cluster ID for each sample
        public double[] wcss; // Within-cluster sum of squares for K=1..10
        public int optimalK;
        public double fstGlobal; // Average Fst between clusters
        public double[][] distanceMatrix; // Pairwise genetic distance
        public int[] dbscanAssignments; // -1 for noise, 0+ for clusters
        public List<double[]> treeSegments; // [x0, y0, x1, y1] for dendrogram
        public int[] gmmAssignments; // GMM cluster labels
        public double[][] dapcMatrix; // Linear Discriminants (LDs)
        public double[][] ancestryProportions; // [numSamples][k] probabilities
        public double[][] kinshipMatrix; // Genomic Relationship Matrix (VanRaden)
    }

    /**
     * Performs PCA on a VCF file.
     * @param vcfFile Path to VCF.
     * @param ploidy Ploidy of the organism.
     * @param numPCs Number of principal components to extract.
     * @param minMaf Minimum MAF for filtering SNPs.
     * @param maxMissing Maximum missingness for filtering SNPs.
     */
    public PcaResult computePCA(String vcfFile, int ploidy, int numPCs, double minMaf, double maxMissing) throws IOException {
        String[] sampleNames = VcfFastReader.getSampleIds(vcfFile);
        int numSamples = sampleNames.length;

        System.out.println("[PCA] Loading and filtering SNPs...");
        List<double[]> dosageMatrix = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
            String line;
            // Identify format index
            int adIdx = -1;
            int gtIdx = -1;

            while ((line = br.readLine()) != null) {
                if (line.startsWith("##FORMAT=<ID=GT")) gtIdx = 0; // default assumption
                if (line.startsWith("#")) continue;

                String[] cols = line.split("\t");
                if (cols.length < 9 + numSamples) continue;

                // Detect indices from FORMAT column
                String[] fmt = cols[8].split(":");
                adIdx = -1; gtIdx = -1;
                for (int i=0; i<fmt.length; i++) {
                    if (fmt[i].equals("GT")) gtIdx = i;
                    if (fmt[i].equals("AD") || fmt[i].equals("BSDP")) adIdx = i;
                }

                double[] siteDosages = new double[numSamples];
                int missingCount = 0;
                double alleleSum = 0;
                int genotypedCount = 0;

                for (int i = 0; i < numSamples; i++) {
                    String gData = cols[9 + i];
                    String[] parts = gData.split(":");
                    double dosage = -1;

                    // Try GT first
                    if (gtIdx != -1 && parts.length > gtIdx && !parts[gtIdx].startsWith(".")) {
                        String[] alleles = parts[gtIdx].split("[/|]");
                        int altCount = 0;
                        for (String a : alleles) if (!a.equals("0")) altCount++;
                        dosage = (double) altCount; // raw dosage (0 to ploidy)
                    } 
                    // Fallback to AD
                    else if (adIdx != -1 && parts.length > adIdx && !parts[adIdx].equals(".")) {
                        String[] ad = parts[adIdx].split(",");
                        if (ad.length >= 2) {
                            try {
                                double r = Double.parseDouble(ad[0]);
                                double a = Double.parseDouble(ad[1]);
                                if (r + a >= 5) dosage = (a / (r + a)) * ploidy;
                            } catch (NumberFormatException ignored) {}
                        }
                    }

                    if (dosage == -1) {
                        missingCount++;
                        siteDosages[i] = Double.NaN;
                    } else {
                        siteDosages[i] = dosage;
                        alleleSum += dosage;
                        genotypedCount++;
                    }
                }

                double missingFreq = (double) missingCount / numSamples;
                if (missingFreq > maxMissing || genotypedCount == 0) continue;

                double freq = alleleSum / (genotypedCount * ploidy);
                double maf = Math.min(freq, 1.0 - freq);
                if (maf < minMaf) continue;

                // Simple mean imputation for missing values
                double meanDosage = alleleSum / genotypedCount;
                for (int i = 0; i < numSamples; i++) {
                    if (Double.isNaN(siteDosages[i])) siteDosages[i] = meanDosage;
                }

                dosageMatrix.add(siteDosages);
            }
        }

        int numMarkers = dosageMatrix.size();
        System.out.println("[PCA] Matrix dimensions: " + numMarkers + " SNPs x " + numSamples + " Samples");

        if (numMarkers < 2) {
            throw new IOException("Not enough SNPs passed filters for PCA.");
        }

        // Build EJML Matrix (Markers as rows, Samples as columns)
        // For genomic PCA, we usually work with M = numSamples x numMarkers
        DMatrixRMaj matrix = new DMatrixRMaj(numSamples, numMarkers);
        for (int j = 0; j < numMarkers; j++) {
            double[] markerData = dosageMatrix.get(j);
            // Calculate mean and variance for standardization
            double sum = 0;
            for (double d : markerData) sum += d;
            double mean = sum / numSamples;
            
            double p = mean / ploidy;
            double std = Math.sqrt(ploidy * p * (1.0 - p));
            if (std == 0) std = 1.0;

            for (int i = 0; i < numSamples; i++) {
                matrix.set(i, j, (markerData[i] - mean) / std);
            }
        }

        System.out.println("[PCA] Computing Singular Value Decomposition (SVD)...");
        // needU=true, needV=true, compact=true
        SingularValueDecomposition_F64<DMatrixRMaj> svd = DecompositionFactory_DDRM.svd(matrix.numRows, matrix.numCols, true, true, true);
        if (!svd.decompose(matrix)) {
            throw new RuntimeException("SVD decomposition failed.");
        }

        DMatrixRMaj V = svd.getV(null, false); // Eigenvectors (Right singular vectors)
        DMatrixRMaj U = svd.getU(null, false); // Loadings (Left singular vectors)
        double[] singularValues = svd.getSingularValues();

        PcaResult result = new PcaResult();
        result.sampleNames = sampleNames;
        int actualPCs = Math.min(numPCs, singularValues.length);
        result.pcMatrix = new double[numSamples][actualPCs];
        result.eigenvalues = new double[actualPCs];
        
        double totalVariance = 0;
        for (double s : singularValues) totalVariance += (s * s);

        result.explainedVariance = new double[actualPCs];
        for (int j = 0; j < actualPCs; j++) {
            result.eigenvalues[j] = singularValues[j] * singularValues[j];
            result.explainedVariance[j] = result.eigenvalues[j] / totalVariance;
            for (int i = 0; i < numSamples; i++) {
                result.pcMatrix[i][j] = U.get(i, j) * singularValues[j];
            }
        }

        // Perform K-Means clustering to find populations
        System.out.println("[PCA] Estimating optimal number of populations (K)...");
        result.wcss = new double[10];
        for (int k = 1; k <= 10; k++) {
            result.wcss[k-1] = runKMeans(result.pcMatrix, k, 50, null);
        }
        
        // Find "elbow" (simple version: find K where reduction slows down)
        result.optimalK = 1;
        for (int k = 1; k < 9; k++) {
            double drop = (result.wcss[k-1] - result.wcss[k]) / result.wcss[0];
            if (drop > 0.1) result.optimalK = k + 1; // Significant drop
        }
        
        System.out.println("[PCA] Optimal K estimated: " + result.optimalK);
        result.clusterAssignments = new int[numSamples];
        runKMeans(result.pcMatrix, result.optimalK, 100, result.clusterAssignments);

        // Compute Global Fst between clusters
        if (result.optimalK > 1) {
            double sumFst = 0;
            int countFst = 0;
            for (double[] dosages : dosageMatrix) {
                double fst = calculateFst(dosages, result.clusterAssignments, result.optimalK, ploidy);
                if (!Double.isNaN(fst)) {
                    sumFst += fst;
                    countFst++;
                }
            }
            result.fstGlobal = countFst > 0 ? sumFst / countFst : 0;
            System.out.println("[PCA] Global Fst between clusters: " + String.format("%.4f", result.fstGlobal));
        }

        // Compute Distance Matrix (Euclidean on PCs)
        result.distanceMatrix = new double[numSamples][numSamples];
        for (int i = 0; i < numSamples; i++) {
            for (int j = i + 1; j < numSamples; j++) {
                double dist = 0;
                for (int l = 0; l < actualPCs; l++) {
                    dist += Math.pow(result.pcMatrix[i][l] - result.pcMatrix[j][l], 2);
                }
                dist = Math.sqrt(dist);
                result.distanceMatrix[i][j] = dist;
                result.distanceMatrix[j][i] = dist;
            }
        }

        // Perform DBSCAN to find noise/outliers
        // Dynamic EPS estimation: use the 10th percentile of distances as a heuristic
        double eps = estimateEps(result.distanceMatrix);
        System.out.println("[PCA] Running DBSCAN (eps=" + String.format("%.2f", eps) + ", minPts=3)...");
        result.dbscanAssignments = runDBSCAN(result.distanceMatrix, eps, 3);

        // Compute Hierarchical Tree (UPGMA)
        System.out.println("[PCA] Computing Hierarchical Tree (UPGMA)...");
        result.treeSegments = computeUPGMA(result.distanceMatrix);

        // Run GMM
        System.out.println("[PCA] Running GMM Clustering (K=" + result.optimalK + ")...");
        result.gmmAssignments = runGMM(result.pcMatrix, result.optimalK, result);

        // Run DAPC
        System.out.println("[PCA] Running DAPC (Discriminant Analysis on PCs)...");
        result.dapcMatrix = runDAPC(result.pcMatrix, result.clusterAssignments, result.optimalK);

        // Compute Kinship Matrix (VanRaden)
        System.out.println("[PCA] Computing Genomic Relationship Matrix (Kinship - VanRaden)...");
        result.kinshipMatrix = computeKinshipVanRaden(dosageMatrix, ploidy);

        return result;
    }

    private double[][] computeKinshipVanRaden(List<double[]> dosageMatrix, int ploidy) {
        int m = dosageMatrix.size();
        int n = dosageMatrix.get(0).length;
        double[][] kinship = new double[n][n];
        
        // 1. Calculate frequencies and scale factor
        double scale = 0;
        double[] p = new double[m];
        for (int j = 0; j < m; j++) {
            double[] dosages = dosageMatrix.get(j);
            double sum = 0;
            int count = 0;
            for (double d : dosages) {
                if (!Double.isNaN(d)) {
                    sum += d;
                    count++;
                }
            }
            p[j] = sum / (count * ploidy);
            scale += ploidy * p[j] * (1.0 - p[j]);
        }

        // 2. Compute Z * Z' where Z is centered genotypes
        // Z_ji = dosage_ji - ploidy * p_j
        // Kinship_i1,i2 = sum_j (Z_ji1 * Z_ji2) / scale
        for (int i1 = 0; i1 < n; i1++) {
            for (int i2 = i1; i2 < n; i2++) {
                double sumZ = 0;
                for (int j = 0; j < m; j++) {
                    double z1 = dosageMatrix.get(j)[i1] - (ploidy * p[j]);
                    double z2 = dosageMatrix.get(j)[i2] - (ploidy * p[j]);
                    sumZ += z1 * z2;
                }
                double val = sumZ / scale;
                kinship[i1][i2] = val;
                kinship[i2][i1] = val;
            }
        }
        
        return kinship;
    }

    private double[][] runDAPC(double[][] pcs, int[] groups, int k) {
        int n = pcs.length;
        int d = Math.min(pcs[0].length, 30); // Use top 30 PCs to avoid overfitting
        
        // 1. Calculate Group Means
        double[][] groupMeans = new double[k][d];
        int[] counts = new int[k];
        for (int i = 0; i < n; i++) {
            int g = groups[i];
            counts[g]++;
            for (int j = 0; j < d; j++) groupMeans[g][j] += pcs[i][j];
        }
        for (int g = 0; g < k; g++) {
            for (int j = 0; j < d; j++) groupMeans[g][j] /= Math.max(1, counts[g]);
        }

        // 2. Center data per group (Within-group scatter)
        double[][] centered = new double[n][d];
        for (int i = 0; i < n; i++) {
            int g = groups[i];
            for (int j = 0; j < d; j++) centered[i][j] = pcs[i][j] - groupMeans[g][j];
        }

        // 3. For DAPC, we often perform PCA on the group means
        // But a more robust way is to project onto the discriminant axes.
        // Simplified DAPC: PCA on Group Means in the PC space
        double[][] meanScatter = new double[k][d];
        for (int g = 0; g < k; g++) {
            for (int j = 0; j < d; j++) meanScatter[g][j] = groupMeans[g][j];
        }

        // Use SVD on group means to find the discriminant axes (LDs)
        double[][] ldAxes = computeSVD(meanScatter, Math.min(k - 1, d));
        
        // Project all samples onto LDs
        double[][] ldProjected = new double[n][ldAxes.length];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < ldAxes.length; j++) {
                double val = 0;
                for (int l = 0; l < d; l++) val += pcs[i][l] * ldAxes[j][l];
                ldProjected[i][j] = val;
            }
        }
        
        return ldProjected;
    }

    private double[][] computeSVD(double[][] matrix, int components) {
        // Simplified SVD via Power Iteration for top eigenvectors
        int d = matrix[0].length;
        int n = matrix.length;
        double[][] eigenvectors = new double[components][d];
        double[][] current = new double[n][d];
        for(int i=0; i<n; i++) System.arraycopy(matrix[i], 0, current[i], 0, d);

        for (int c = 0; c < components; c++) {
            double[] v = new double[d];
            Random rnd = new Random(42 + c);
            for (int i = 0; i < d; i++) v[i] = rnd.nextDouble();
            
            for (int iter = 0; iter < 50; iter++) {
                double[] nextV = new double[d];
                for (int i = 0; i < n; i++) {
                    double dot = 0;
                    for (int j = 0; j < d; j++) dot += current[i][j] * v[j];
                    for (int j = 0; j < d; j++) nextV[j] += current[i][j] * dot;
                }
                double norm = 0;
                for (double val : nextV) norm += val * val;
                norm = Math.sqrt(norm);
                if (norm < 1e-10) break;
                for (int i = 0; i < d; i++) v[i] = nextV[i] / norm;
            }
            eigenvectors[c] = v;
            
            // Deflation
            for (int i = 0; i < n; i++) {
                double dot = 0;
                for (int j = 0; j < d; j++) dot += current[i][j] * v[j];
                for (int j = 0; j < d; j++) current[i][j] -= dot * v[j];
            }
        }
        return eigenvectors;
    }

    private int[] runGMM(double[][] data, int k, PcaResult result) {
        int n = data.length;
        int d = Math.min(data[0].length, 3);
        
        double[][] means = new double[k][d];
        double[] weights = new double[k];
        double[][] variances = new double[k][d];
        
        Random rnd = new Random(42);
        for (int i = 0; i < k; i++) {
            means[i] = data[rnd.nextInt(n)].clone();
            weights[i] = 1.0 / k;
            for (int j = 0; j < d; j++) variances[i][j] = 1.0;
        }

        double[][] resp = new double[n][k];
        for (int iter = 0; iter < 50; iter++) {
            // E-Step
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < k; j++) {
                    double prob = weights[j];
                    for (int dim = 0; dim < d; dim++) {
                        double diff = data[i][dim] - means[j][dim];
                        prob *= Math.exp(-0.5 * diff * diff / (variances[j][dim] + 1e-6)) / (Math.sqrt(2 * Math.PI * variances[j][dim]) + 1e-6);
                    }
                    resp[i][j] = prob;
                    sum += prob;
                }
                if (sum > 0) {
                    for (int j = 0; j < k; j++) resp[i][j] /= sum;
                }
            }

            // M-Step
            for (int j = 0; j < k; j++) {
                double nK = 0;
                for (int i = 0; i < n; i++) nK += resp[i][j];
                
                if (nK > 1e-4) {
                    weights[j] = nK / n;
                    for (int dim = 0; dim < d; dim++) {
                        double m = 0;
                        for (int i = 0; i < n; i++) m += resp[i][j] * data[i][dim];
                        means[j][dim] = m / nK;
                        
                        double v = 0;
                        for (int i = 0; i < n; i++) {
                            double diff = data[i][dim] - means[j][dim];
                            v += resp[i][j] * diff * diff;
                        }
                        variances[j][dim] = Math.max(v / nK, 1e-4);
                    }
                }
            }
        }

        result.ancestryProportions = resp; // Save for barplot
        
        int[] assignments = new int[n];
        for (int i = 0; i < n; i++) {
            int best = 0;
            for (int j = 1; j < k; j++) {
                if (resp[i][j] > resp[i][best]) best = j;
            }
            assignments[i] = best;
        }
        return assignments;
    }

    private List<double[]> computeUPGMA(double[][] dists) {
        int n = dists.length;
        List<double[]> segments = new ArrayList<>();
        
        // Node representation
        class Node {
            int id;
            double x; // Leaf position (0 to n-1)
            double height;
            int size;
            Node left, right;
            Node(int id, double x) { this.id = id; this.x = x; this.size = 1; this.height = 0; }
        }

        List<Node> activeNodes = new ArrayList<>();
        for (int i = 0; i < n; i++) activeNodes.add(new Node(i, i));

        double[][] currentDists = new double[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(dists[i], 0, currentDists[i], 0, n);

        while (activeNodes.size() > 1) {
            // Find min dist
            int bestI = 0, bestJ = 1;
            double minDist = Double.MAX_VALUE;
            for (int i = 0; i < activeNodes.size(); i++) {
                for (int j = i + 1; j < activeNodes.size(); j++) {
                    if (currentDists[i][j] < minDist) {
                        minDist = currentDists[i][j];
                        bestI = i; bestJ = j;
                    }
                }
            }

            Node n1 = activeNodes.get(bestI);
            Node n2 = activeNodes.get(bestJ);
            double mergeHeight = minDist / 2.0;

            // New node in the middle
            Node newNode = new Node(-1, (n1.x + n2.x) / 2.0);
            newNode.height = mergeHeight;
            newNode.size = n1.size + n2.size;
            newNode.left = n1;
            newNode.right = n2;

            // Add segments for Plotly [x0, y0, x1, y1]
            // Horizontal bar
            segments.add(new double[]{n1.x, mergeHeight, n2.x, mergeHeight});
            // Vertical bars to children
            segments.add(new double[]{n1.x, n1.height, n1.x, mergeHeight});
            segments.add(new double[]{n2.x, n2.height, n2.x, mergeHeight});

            // Update distance matrix (UPGMA: average distance)
            double[] newNodeDists = new double[activeNodes.size()];
            for (int k = 0; k < activeNodes.size(); k++) {
                if (k == bestI || k == bestJ) continue;
                // UPGMA formula: d(k, i+j) = (d(k,i)*size(i) + d(k,j)*size(j)) / (size(i)+size(j))
                newNodeDists[k] = (currentDists[bestI][k] * n1.size + currentDists[bestJ][k] * n2.size) / (double)(n1.size + n2.size);
            }

            // Remove n1, n2 and add newNode
            // We'll replace bestI with newNode and remove bestJ
            activeNodes.set(bestI, newNode);
            activeNodes.remove(bestJ);

            // Update currentDists matrix
            double[][] nextDists = new double[activeNodes.size()][activeNodes.size()];
            for (int i = 0, oldI = 0; i < nextDists.length; i++, oldI++) {
                if (oldI == bestJ) oldI++; // Skip removed index
                for (int j = 0, oldJ = 0; j < nextDists.length; j++, oldJ++) {
                    if (oldJ == bestJ) oldJ++; // Skip removed index
                    
                    if (i == bestI) nextDists[i][j] = newNodeDists[oldJ];
                    else if (j == bestI) nextDists[i][j] = newNodeDists[oldI];
                    else nextDists[i][j] = currentDists[oldI][oldJ];
                }
            }
            currentDists = nextDists;
        }

        return segments;
    }

    private double estimateEps(double[][] dists) {
        List<Double> allDists = new ArrayList<>();
        for (int i = 0; i < dists.length; i++) {
            for (int j = i + 1; j < dists.length; j++) {
                allDists.add(dists[i][j]);
            }
        }
        Collections.sort(allDists);
        // Use the 5th percentile distance as epsilon for a tight but inclusive clustering
        int index = (int) (allDists.size() * 0.05);
        return allDists.get(index);
    }

    private int[] runDBSCAN(double[][] dists, double eps, int minPts) {
        int n = dists.length;
        int[] assignments = new int[n];
        Arrays.fill(assignments, -2); // -2: unvisited
        int clusterId = 0;

        for (int i = 0; i < n; i++) {
            if (assignments[i] != -2) continue;

            List<Integer> neighbors = getNeighbors(dists, i, eps);
            if (neighbors.size() < minPts) {
                assignments[i] = -1; // Noise
            } else {
                assignments[i] = clusterId;
                expandCluster(dists, assignments, neighbors, clusterId, eps, minPts);
                clusterId++;
            }
        }
        return assignments;
    }

    private void expandCluster(double[][] dists, int[] assignments, List<Integer> neighbors, int clusterId, double eps, int minPts) {
        for (int i = 0; i < neighbors.size(); i++) {
            int p = neighbors.get(i);
            if (assignments[p] == -1) assignments[p] = clusterId;
            if (assignments[p] != -2) continue;

            assignments[p] = clusterId;
            List<Integer> nextNeighbors = getNeighbors(dists, p, eps);
            if (nextNeighbors.size() >= minPts) {
                for (int next : nextNeighbors) {
                    if (!neighbors.contains(next)) neighbors.add(next);
                }
            }
        }
    }

    private List<Integer> getNeighbors(double[][] dists, int i, double eps) {
        List<Integer> neighbors = new ArrayList<>();
        for (int j = 0; j < dists.length; j++) {
            if (dists[i][j] <= eps) neighbors.add(j);
        }
        return neighbors;
    }

    private double calculateFst(double[] dosages, int[] assignments, int k, int ploidy) {
        double[] p = new double[k];
        int[] n = new int[k];
        double totalSum = 0;
        int totalN = 0;

        for (int i = 0; i < dosages.length; i++) {
            if (Double.isNaN(dosages[i])) continue;
            p[assignments[i]] += dosages[i];
            n[assignments[i]]++;
            totalSum += dosages[i];
            totalN++;
        }

        if (totalN == 0) return Double.NaN;
        
        double pTotal = totalSum / (totalN * ploidy);
        double ht = 2 * pTotal * (1.0 - pTotal);
        if (ht < 0.0001) return 0;

        double hs = 0;
        for (int i = 0; i < k; i++) {
            if (n[i] > 0) {
                double pi = p[i] / (n[i] * ploidy);
                double hi = 2 * pi * (1.0 - pi);
                hs += ((double) n[i] / totalN) * hi;
            }
        }

        return (ht - hs) / ht;
    }

    private double runKMeans(double[][] data, int k, int maxIter, int[] assignments) {
        int n = data.length;
        int d = data[0].length;
        double[][] centroids = new double[k][d];
        int[] localAssignments = new int[n];
        
        // Init centroids (random samples)
        Random rnd = new Random(42);
        for (int i = 0; i < k; i++) {
            System.arraycopy(data[rnd.nextInt(n)], 0, centroids[i], 0, d);
        }

        double totalDist = 0;
        for (int iter = 0; iter < maxIter; iter++) {
            totalDist = 0;
            // E-step: Assign
            for (int i = 0; i < n; i++) {
                double minDist = Double.MAX_VALUE;
                int bestK = 0;
                for (int c = 0; c < k; c++) {
                    double dist = 0;
                    for (int l = 0; l < d; l++) dist += Math.pow(data[i][l] - centroids[c][l], 2);
                    if (dist < minDist) { minDist = dist; bestK = c; }
                }
                localAssignments[i] = bestK;
                totalDist += minDist;
            }
            // M-step: Update centroids
            double[][] nextCentroids = new double[k][d];
            int[] counts = new int[k];
            for (int i = 0; i < n; i++) {
                int c = localAssignments[i];
                counts[c]++;
                for (int l = 0; l < d; l++) nextCentroids[c][l] += data[i][l];
            }
            for (int c = 0; c < k; c++) {
                if (counts[c] > 0) {
                    for (int l = 0; l < d; l++) centroids[c][l] = nextCentroids[c][l] / counts[c];
                }
            }
        }
        
        if (assignments != null) System.arraycopy(localAssignments, 0, assignments, 0, n);
        return totalDist;
    }

    public void exportPca(PcaResult result, String outputPath) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputPath))) {
            // Header
            pw.print("Sample,KMeans_Cluster,DBSCAN_Cluster,GMM_Cluster");
            for (int j = 0; j < result.pcMatrix[0].length; j++) pw.print(",PC" + (j + 1));
            for (int j = 0; j < result.dapcMatrix[0].length; j++) pw.print(",LD" + (j + 1));
            for (int j = 0; j < result.ancestryProportions[0].length; j++) pw.print(",Ancestry_Q" + (j + 1));
            pw.println();

            // Data
            for (int i = 0; i < result.sampleNames.length; i++) {
                String dbscanStr = result.dbscanAssignments[i] == -1 ? "Noise" : String.valueOf(result.dbscanAssignments[i] + 1);
                pw.print(result.sampleNames[i] + "," + (result.clusterAssignments[i] + 1) + "," + dbscanStr + "," + (result.gmmAssignments[i] + 1));
                
                // PCA
                for (int j = 0; j < result.pcMatrix[0].length; j++) pw.printf(Locale.US, ",%.6f", result.pcMatrix[i][j]);
                // DAPC
                for (int j = 0; j < result.dapcMatrix[0].length; j++) pw.printf(Locale.US, ",%.6f", result.dapcMatrix[i][j]);
                // Ancestry
                for (int j = 0; j < result.ancestryProportions[0].length; j++) pw.printf(Locale.US, ",%.6f", result.ancestryProportions[i][j]);
                
                pw.println();
            }
        }
        System.out.println("[PCA] Results exported to: " + outputPath);
    }

    public void exportKinship(PcaResult result, String outputPath) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputPath))) {
            // Header: Sample names
            pw.print("Sample");
            for (String name : result.sampleNames) pw.print("," + name);
            pw.println();

            // Rows
            for (int i = 0; i < result.sampleNames.length; i++) {
                pw.print(result.sampleNames[i]);
                for (int j = 0; j < result.sampleNames.length; j++) {
                    pw.printf(Locale.US, ",%.6f", result.kinshipMatrix[i][j]);
                }
                pw.println();
            }
        }
        System.out.println("[PCA] Kinship matrix exported to: " + outputPath);
    }
}
