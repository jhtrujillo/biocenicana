package org.cenicana.bio.io;

import org.cenicana.bio.core.GeneFeature;
import org.cenicana.bio.utils.ChromosomeNormalizer;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Specialized parser for GFF3 files to extract Gene features for the Annotation Dashboard.
 */
public class GeneAnnotationParser {

    public static Map<String, List<GeneFeature>> parseGenes(String gffPath) throws IOException {
        Map<String, List<GeneFeature>> geneIndex = new HashMap<>();
        Map<String, GeneFeature> genesById = new HashMap<>();
        Map<String, String> parentMap = new HashMap<>(); // ID -> ParentID

        try (BufferedReader br = new BufferedReader(new FileReader(gffPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                
                String[] cols = line.split("\t");
                if (cols.length < 9) continue;
                
                String type = cols[2].toLowerCase();
                String chr = ChromosomeNormalizer.normalize(cols[0]);
                long start = Long.parseLong(cols[3]);
                long end = Long.parseLong(cols[4]);
                String strand = cols[6];
                String attributes = cols[8];

                String id = getAttribute(attributes, "ID");
                String parent = getAttribute(attributes, "Parent");

                if (type.equals("gene") || ((type.equals("mrna") || type.equals("transcript")) && parent == null)) {
                    GeneFeature gene = new GeneFeature(chr, start, end, strand);
                    parseAttributes(gene, attributes);
                    if (id != null) genesById.put(id, gene);
                    geneIndex.computeIfAbsent(chr, k -> new ArrayList<>()).add(gene);
                } else if (type.equals("mrna") || type.equals("transcript")) {
                    if (id != null && parent != null) {
                        parentMap.put(id, parent);
                        GeneFeature parentGene = genesById.get(parent);
                        if (parentGene != null) parseAttributes(parentGene, attributes);
                    }
                } else if (type.equals("exon") || type.equals("cds") || type.equals("five_prime_utr") || type.equals("three_prime_utr")) {
                    if (parent != null) {
                        // Find the ultimate gene parent
                        String geneId = parent;
                        if (parentMap.containsKey(geneId)) geneId = parentMap.get(geneId);
                        
                        GeneFeature gene = genesById.get(geneId);
                        if (gene != null) {
                            gene.addSubFeature(type, start, end);
                        }
                    }
                }
            }
        }
        
        for (List<GeneFeature> genes : geneIndex.values()) {
            genes.sort(Comparator.comparingLong(GeneFeature::getStart));
        }
        
        return geneIndex;
    }

    private static String getAttribute(String attrStr, String key) {
        for (String pair : attrStr.split(";")) {
            String[] kv = pair.split("=");
            if (kv.length == 2 && kv[0].equalsIgnoreCase(key)) return kv[1];
        }
        return null;
    }

    private static void parseAttributes(GeneFeature gene, String attrStr) {
        String[] pairs = attrStr.split(";");
        for (String pair : pairs) {
            String[] kv = pair.split("=");
            if (kv.length != 2) continue;
            
            String key = kv[0].trim();
            String value = kv[1].trim();
            
            switch (key.toLowerCase()) {
                case "id": gene.setId(value); break;
                case "name": gene.setName(value); break;
                case "note":
                case "description":
                case "product":
                    gene.setNote(value.replace("%20", " ").replace("+", " "));
                    break;
                case "ontology_term":
                case "go":
                    for (String v : value.split(",")) gene.addGoTerm(v);
                    break;
                case "dbxref":
                    parseDbxref(gene, value);
                    break;
                default:
                    gene.getOtherAttributes().put(key, value);
                    break;
            }
        }
    }

    private static void parseDbxref(GeneFeature gene, String value) {
        for (String part : value.split(",")) {
            String[] sub = part.split(":");
            if (sub.length != 2) continue;
            String db = sub[1].toLowerCase(); // Often format is "Database:ID"
            String val = sub[1];
            
            // Re-evaluating the split logic for Dbxref
            String dbName = sub[0].toLowerCase();
            if (dbName.equals("pfam")) gene.addPfamDomain(val);
            else if (dbName.equals("interpro")) gene.addInterproDomain(val);
            else if (dbName.equals("go")) gene.addGoTerm(val);
        }
    }
}
