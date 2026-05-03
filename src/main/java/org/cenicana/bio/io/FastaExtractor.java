package org.cenicana.bio.io;

import java.io.*;
import java.util.*;

/**
 * Utility to quickly extract sequences from large FASTA files using an index.
 */
public class FastaExtractor {
    private File fastaFile;
    private Map<String, Long> index = new HashMap<>();

    public FastaExtractor(File file) throws IOException {
        this.fastaFile = file;
        buildIndex();
    }

    private void buildIndex() throws IOException {
        try (RandomAccessFile raf = new RandomAccessFile(fastaFile, "r")) {
            long pos = 0;
            String line;
            while ((line = raf.readLine()) != null) {
                if (line.startsWith(">")) {
                    // Extract ID and normalize it
                    String rawId = line.substring(1).split("\\s+")[0];
                    String normalizedId = org.cenicana.bio.utils.ChromosomeNormalizer.normalize(rawId);
                    index.put(normalizedId, pos);
                }
                pos = raf.getFilePointer();
            }
        }
    }

    public String extract(String id) {
        return extract(id, 1, -1);
    }

    public String extract(String id, long start, long end) {
        Long pos = index.get(id);
        if (pos == null) return null;

        StringBuilder seq = new StringBuilder();
        try (RandomAccessFile raf = new RandomAccessFile(fastaFile, "r")) {
            raf.seek(pos);
            raf.readLine(); // skip header
            
            // For sub-sequences, we could optimize more, but this is robust
            String line;
            long currentPos = 1;
            while ((line = raf.readLine()) != null) {
                if (line.startsWith(">")) break;
                line = line.trim();
                int len = line.length();
                
                if (end == -1 || (currentPos + len >= start && currentPos <= end)) {
                    for (int i = 0; i < len; i++) {
                        long globalPos = currentPos + i;
                        if (globalPos >= start && (end == -1 || globalPos <= end)) {
                            seq.append(line.charAt(i));
                        }
                        if (end != -1 && globalPos > end) break;
                    }
                }
                currentPos += len;
                if (end != -1 && currentPos > end) break;
            }
        } catch (IOException e) {
            return null;
        }
        return seq.toString();
    }
}
