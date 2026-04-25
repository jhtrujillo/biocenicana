package org.cenicana.bio.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Reader for FASTA files (CDS, Protein, etc.)
 */
public class FastaReader {

    /**
     * Reads a FASTA file and returns a map of ID -> Sequence.
     */
    public Map<String, String> read(String filePath) throws IOException {
        Map<String, String> sequences = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            String currentId = null;
            StringBuilder currentSeq = new StringBuilder();

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    if (currentId != null) {
                        sequences.put(currentId, currentSeq.toString());
                    }
                    // Extract ID before first space
                    currentId = line.substring(1).split("\\s+")[0];
                    currentSeq = new StringBuilder();
                } else {
                    currentSeq.append(line);
                }
            }
            if (currentId != null) {
                sequences.put(currentId, currentSeq.toString());
            }
        }
        return sequences;
    }
}
