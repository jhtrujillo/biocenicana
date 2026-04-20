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
                if (line.startsWith("#CHROM")) {
                    String[] header = line.split("\t");
                    if (header.length > 8) {
                        String[] fmt = header[8].split(":"); // This is wrong, header[8] is "FORMAT" literal
                    }
                    continue;
                }
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

        return result;
    }

    public void exportPca(PcaResult result, String outputPath) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputPath))) {
            pw.print("Sample");
            for (int j = 0; j < result.pcMatrix[0].length; j++) pw.print(",PC" + (j + 1));
            pw.println();
            for (int i = 0; i < result.sampleNames.length; i++) {
                pw.print(result.sampleNames[i]);
                for (int j = 0; j < result.pcMatrix[0].length; j++) {
                    pw.printf(Locale.US, ",%.6f", result.pcMatrix[i][j]);
                }
                pw.println();
            }
        }
        System.out.println("[PCA] Results exported to: " + outputPath);
    }
}
