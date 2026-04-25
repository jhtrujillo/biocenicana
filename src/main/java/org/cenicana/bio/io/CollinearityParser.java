package org.cenicana.bio.io;

import org.cenicana.bio.model.SyntenicBlock;
import org.cenicana.bio.model.SyntenicPair;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Parser for McScanX .collinearity files.
 */
public class CollinearityParser {

    public List<SyntenicBlock> parse(String filePath) throws IOException {
        List<SyntenicBlock> blocks = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            SyntenicBlock currentBlock = null;
            boolean isMcScanX = false;
            boolean checkedFormat = false;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (!checkedFormat) {
                    if (line.startsWith("##") || line.contains("Alignment")) {
                        isMcScanX = true;
                    }
                    checkedFormat = true;
                }

                if (isMcScanX) {
                    if (line.startsWith("## Alignment")) {
                        String[] parts = line.split("\\s+");
                        String blockId = parts[2].replace(":", "");
                        currentBlock = new SyntenicBlock(blockId);
                        for (String p : parts) {
                            if (p.startsWith("score=")) currentBlock.setScore(Double.parseDouble(p.split("=")[1]));
                            else if (p.startsWith("e_value=")) currentBlock.seteValue(Double.parseDouble(p.split("=")[1]));
                            else if (p.startsWith("N=")) currentBlock.setNumGenes(Integer.parseInt(p.split("=")[1]));
                            else if (p.contains("plus")) currentBlock.setOrientation("plus");
                            else if (p.contains("minus")) currentBlock.setOrientation("minus");
                        }
                        blocks.add(currentBlock);
                    } else if (line.contains(":") && !line.startsWith("#")) {
                        String[] mainParts = line.split(":");
                        if (mainParts.length >= 2) {
                            String[] dataParts = mainParts[1].trim().split("\\s+");
                            if (dataParts.length >= 3 && currentBlock != null) {
                                String g1 = dataParts[0];
                                String g2 = dataParts[1];
                                double ev = parseDoubleSafe(dataParts[2]);
                                currentBlock.addPair(new SyntenicPair(g1, g2, ev));
                            }
                        }
                    }
                } else {
                    // Generic TSV / SynMap2 format
                    // Expected: BlockID \t Gene1 \t Gene2 \t Score \t EValue
                    if (line.startsWith("#")) continue;
                    String[] parts = line.split("\t");
                    if (parts.length >= 3) {
                        String bId = parts[0];
                        String g1 = parts[1];
                        String g2 = parts[2];
                        double ev = parts.length > 4 ? parseDoubleSafe(parts[4]) : 0;
                        
                        // Check if we need to start a new block
                        if (currentBlock == null || !currentBlock.getBlockId().equals(bId)) {
                            currentBlock = new SyntenicBlock(bId);
                            blocks.add(currentBlock);
                        }
                        currentBlock.addPair(new SyntenicPair(g1, g2, ev));
                    }
                }
            }
        }
        return blocks;
    }

    private double parseDoubleSafe(String s) {
        try {
            return Double.parseDouble(s);
        } catch (Exception e) {
            return 0;
        }
    }
}
