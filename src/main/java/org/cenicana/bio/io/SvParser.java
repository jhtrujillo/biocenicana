package org.cenicana.bio.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Parser for Structural Variants (SVs) from VCF files.
 * Identifies regions with large deletions, insertions, or inversions.
 */
public class SvParser {

    public static class SvRegion {
        public String chromosome;
        public long start;
        public long end;
        public String type;

        public SvRegion(String chr, long start, long end, String type) {
            this.chromosome = chr;
            this.start = start;
            this.end = end;
            this.type = type;
        }
    }

    public List<SvRegion> parse(String vcfPath) throws IOException {
        List<SvRegion> svs = new ArrayList<>();
        if (vcfPath == null || vcfPath.isBlank()) return svs;

        try (BufferedReader br = Files.newBufferedReader(Paths.get(vcfPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.startsWith("#")) continue;

                String[] cols = line.split("\t");
                if (cols.length < 8) continue;

                String chr = cols[0];
                long start = Long.parseLong(cols[1]);
                String info = cols[7];

                // Detect SVTYPE and END
                String type = "SV";
                long end = start;

                if (info.contains("SVTYPE=")) {
                    type = getValue(info, "SVTYPE");
                }
                
                if (info.contains("END=")) {
                    try {
                        end = Long.parseLong(getValue(info, "END"));
                    } catch (NumberFormatException e) {
                        // Fallback: if no END, use length of REF/ALT if it's a simple indel
                        end = start + Math.max(cols[3].length(), cols[4].length());
                    }
                } else {
                    // Estimate end for non-symbolic alleles
                    end = start + Math.max(cols[3].length(), cols[4].length());
                }

                // Only consider SVs (size > 50bp as a heuristic)
                if (Math.abs(end - start) > 50 || !"SV".equals(type)) {
                    svs.add(new SvRegion(chr, start, end, type));
                }
            }
        }
        return svs;
    }

    private String getValue(String info, String key) {
        int idx = info.indexOf(key + "=");
        if (idx == -1) return "NA";
        int start = idx + key.length() + 1;
        int end = info.indexOf(";", start);
        if (end == -1) end = info.length();
        return info.substring(start, end);
    }
}
