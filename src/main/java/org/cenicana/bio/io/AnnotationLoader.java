package org.cenicana.bio.io;

import java.io.*;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Loads functional annotations from either a GFF3 file (column 9 attributes)
 * or a TSV file (GeneID \t Description [optional: \t GO_terms]).
 *
 * Returns a Map<GeneID, Description>.
 */
public class AnnotationLoader {

    public Map<String, String> load(String filePath) throws IOException {
        if (filePath == null || filePath.isBlank()) return new HashMap<>();
        if (filePath.toLowerCase().endsWith(".gff") || filePath.toLowerCase().endsWith(".gff3")) {
            return loadFromGff3(filePath);
        }
        return loadFromTsv(filePath);
    }

    /** Loads GO terms: returns Map<GeneID, List<GO_IDs>> */
    public Map<String, List<String>> loadGoTerms(String filePath) throws IOException {
        if (filePath == null || filePath.isBlank()) return new HashMap<>();
        Map<String, List<String>> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) continue;
                String[] p = line.split("\t");
                
                String id = null;
                String rawGo = null;

                if (p.length >= 9 && (filePath.toLowerCase().endsWith(".gff") || filePath.toLowerCase().endsWith(".gff3"))) {
                    // GFF3 mode
                    String type = p[2].toLowerCase();
                    if (!type.equals("gene") && !type.equals("mrna")) continue;
                    Map<String, String> attrs = parseAttrs(p[8]);
                    if (type.equals("mrna")) {
                        // For mRNA: use Parent gene ID as key (strip isoform suffix)
                        id = attrs.get("Parent");
                        if (id != null) id = id.replaceFirst("^(gene:)", "");
                    }
                    if (id == null) {
                        id = attrs.getOrDefault("ID", attrs.getOrDefault("Name", null));
                        if (id != null) id = id.replaceFirst("^(gene:|mRNA:|transcript:)", "");
                    }
                    if (id != null) {
                        rawGo = firstNonEmpty(attrs, "Ontology_term", "GO", "go", "Dbxref");
                    }
                } else if (p.length >= 3) {
                    // TSV mode: ID \t Desc \t GO1,GO2...
                    id = p[0].trim();
                    rawGo = p[2].trim();
                }

                if (id != null && rawGo != null) {
                    List<String> terms = extractGoIds(rawGo);
                    if (!terms.isEmpty()) {
                        map.computeIfAbsent(id, k -> new ArrayList<>()).addAll(terms);
                    }
                }
            }
        }
        return map;
    }

    private List<String> extractGoIds(String s) {
        List<String> list = new ArrayList<>();
        // Match GO:XXXXXXX patterns
        java.util.regex.Matcher m = java.util.regex.Pattern.compile("GO:\\d+").matcher(s);
        while (m.find()) {
            list.add(m.group());
        }
        return list;
    }

    /** Parse GFF3: extract ID and description/Note from column 9. */
    private Map<String, String> loadFromGff3(String path) throws IOException {
        Map<String, String> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) continue;
                String[] p = line.split("\t");
                if (p.length < 9) continue;
                String type = p[2].toLowerCase();
                if (!type.equals("gene") && !type.equals("mrna")) continue;
                Map<String, String> attrs = parseAttrs(p[8]);
                String id;
                if (type.equals("mrna")) {
                    // Use Parent gene ID so annotation key matches the gene loader
                    id = attrs.get("Parent");
                    if (id != null) id = id.replaceFirst("^(gene:)", "");
                } else {
                    id = attrs.getOrDefault("ID", attrs.getOrDefault("Name", null));
                    if (id != null) id = id.replaceFirst("^(gene:|mRNA:|transcript:)", "");
                }
                if (id == null) continue;
                String desc = firstNonEmpty(attrs,
                        "description", "Description", "Note", "note", "product", "function");
                if (desc != null && !desc.isBlank()) {
                    map.put(id, URLDecoder.decode(desc, StandardCharsets.UTF_8));
                }
            }
        }
        return map;
    }

    /** Parse TSV: GeneID \t Description [\t GO_terms] */
    private Map<String, String> loadFromTsv(String path) throws IOException {
        Map<String, String> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) continue;
                String[] p = line.split("\t");
                if (p.length < 2) continue;
                String id = p[0].trim();
                String desc = p[1].trim();
                if (!id.isBlank() && !desc.isBlank()) map.put(id, desc);
            }
        }
        return map;
    }

    private Map<String, String> parseAttrs(String s) {
        Map<String, String> m = new HashMap<>();
        for (String pair : s.split(";")) {
            String[] kv = pair.split("=", 2);
            if (kv.length == 2) m.put(kv[0].trim(), kv[1].trim());
        }
        return m;
    }

    private String firstNonEmpty(Map<String, String> m, String... keys) {
        for (String k : keys) {
            String v = m.get(k);
            if (v != null && !v.isBlank()) return v;
        }
        return null;
    }
}
