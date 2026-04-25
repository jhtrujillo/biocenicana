package org.cenicana.bio.io;

import org.cenicana.bio.model.Gene;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Parser for GFF3 files.
 */
public class GffParser {

    /**
     * Parses a GFF3 file and returns a list of Gene features.
     * Only lines of type 'gene' or 'mRNA' are typically processed for synteny.
     */
    public List<Gene> parse(String filePath) throws IOException {
        List<Gene> genes = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;

                String[] parts = line.split("\t");
                if (parts.length < 9) continue;

                String chromosome = parts[0];
                String type = parts[2];
                long start = Long.parseLong(parts[3]);
                long end = Long.parseLong(parts[4]);
                String strand = parts[6];
                String attributesStr = parts[8];

                Map<String, String> attributes = parseAttributes(attributesStr);
                String id = attributes.getOrDefault("ID", attributes.getOrDefault("Name", "unknown_" + genes.size()));

                Gene gene = new Gene(id, chromosome, start, end, strand, type);
                for (Map.Entry<String, String> entry : attributes.entrySet()) {
                    gene.addAttribute(entry.getKey(), entry.getValue());
                }
                genes.add(gene);
            }
        }
        return genes;
    }

    private Map<String, String> parseAttributes(String attrStr) {
        Map<String, String> map = new HashMap<>();
        String[] pairs = attrStr.split(";");
        for (String pair : pairs) {
            String[] kv = pair.split("=");
            if (kv.length == 2) {
                map.put(kv[0].trim(), kv[1].trim());
            }
        }
        return map;
    }
}
