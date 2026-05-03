package org.cenicana.bio.io;

import java.io.*;
import java.util.*;

/**
 * Handles reading phenotype files (CSV/TSV) for GWAS analysis.
 */
public class PhenotypeData {
    private Map<String, Double> values = new LinkedHashMap<>();
    private Map<String, Map<String, String>> fixedEffects = new HashMap<>();
    private String traitName;

    /**
     * Loads phenotype data from a file.
     * @param path Path to the CSV or TSV file.
     * @param targetTrait Name of the column containing the trait values.
     * @throws IOException If file cannot be read or trait is not found.
     */
    public void load(String path, String targetTrait, List<String> extraCols) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String header = br.readLine();
            if (header == null) return;
            
            String sep = header.contains(",") ? "," : "\t";
            String[] cols = header.split(sep);
            int idIdx = -1;
            int traitIdx = -1;
            Map<String, Integer> fixedIndices = new HashMap<>();

            for (int i = 0; i < cols.length; i++) {
                String c = cols[i].trim();
                if (c.equalsIgnoreCase("id") || c.equalsIgnoreCase("sample") || c.equalsIgnoreCase("taxa")) {
                    idIdx = i;
                }
                if (c.equalsIgnoreCase(targetTrait)) {
                    traitIdx = i;
                }
                if (extraCols != null) {
                    for (String ex : extraCols) {
                        if (c.equalsIgnoreCase(ex.trim())) fixedIndices.put(ex.trim(), i);
                    }
                }
            }

            if (idIdx == -1) idIdx = 0;
            if (traitIdx == -1) throw new IOException("Trait '" + targetTrait + "' not found.");
            
            this.traitName = cols[traitIdx].trim();

            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(sep);
                if (parts.length > Math.max(idIdx, traitIdx)) {
                    String sampleId = parts[idIdx].trim();
                    String valStr = parts[traitIdx].trim();
                    if (!valStr.isEmpty() && !valStr.equals("NA")) {
                        try {
                            values.put(sampleId, Double.parseDouble(valStr));
                            Map<String, String> fMap = new HashMap<>();
                            for (Map.Entry<String, Integer> entry : fixedIndices.entrySet()) {
                                if (parts.length > entry.getValue()) fMap.put(entry.getKey(), parts[entry.getValue()].trim());
                            }
                            fixedEffects.put(sampleId, fMap);
                        } catch (Exception e) {}
                    }
                }
            }
        }
        System.out.println("[IO] Loaded " + values.size() + " phenotypes for trait: " + traitName);
    }

    public Map<String, String> getFixedEffects(String sampleId) {
        return fixedEffects.getOrDefault(sampleId, new HashMap<>());
    }

    public Double getValue(String sampleId) {
        return values.get(sampleId);
    }

    public String getTraitName() {
        return traitName;
    }

    public Set<String> getSamples() {
        return values.keySet();
    }
    
    /**
     * Returns the values for a specific list of samples, in the same order.
     * Missing samples will have Double.NaN.
     */
    public double[] getOrderedValues(String[] sampleIds) {
        double[] result = new double[sampleIds.length];
        for (int i = 0; i < sampleIds.length; i++) {
            Double v = values.get(sampleIds[i]);
            result[i] = (v != null) ? v : Double.NaN;
        }
        return result;
    }
}
