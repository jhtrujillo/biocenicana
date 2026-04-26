package org.cenicana.bio.cli;

import org.cenicana.bio.core.KaKsCalculator;
import org.cenicana.bio.io.CollinearityParser;
import org.cenicana.bio.io.FastaReader;
import org.cenicana.bio.model.SyntenicBlock;
import org.cenicana.bio.model.SyntenicPair;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.Callable;

@Command(name = "kaks-calc", 
         mixinStandardHelpOptions = true, 
         description = "Calculate Ka/Ks ratios for orthologous gene pairs using Nei-Gojobori method.")
public class KaKsCalcCommand implements Callable<Integer> {

    @Option(names = {"--collinearity"}, required = true, description = "Collinearity file (McScanX output).")
    private String collinearity;

    @Option(names = {"--cds1"}, required = true, description = "CDS FASTA for Genome 1.")
    private String cds1;

    @Option(names = {"--cds2"}, required = true, description = "CDS FASTA for Genome 2.")
    private String cds2;

    @Option(names = {"-o", "--output"}, description = "Output TSV file path.", defaultValue = "kaks_results.tsv")
    private String output;

    @Override
    public Integer call() throws Exception {
        System.out.println("[Phase 1/3] Loading CDS sequences...");
        FastaReader fastaReader = new FastaReader();
        Map<String, String> seqs1 = fastaReader.read(cds1);
        Map<String, String> seqs2 = fastaReader.read(cds2);
        System.out.println("  - Loaded " + seqs1.size() + " sequences from Gen1.");
        System.out.println("  - Loaded " + seqs2.size() + " sequences from Gen2.");

        System.out.println("[Phase 2/3] Parsing collinearity and calculating Ka/Ks...");
        CollinearityParser colParser = new CollinearityParser();
        List<SyntenicBlock> blocks = colParser.parse(collinearity);
        
        KaKsCalculator calculator = new KaKsCalculator();
        int count = 0;

        try (PrintWriter pw = new PrintWriter(new FileWriter(output))) {
            pw.println("#Gene1\tGene2\tKa\tKs\tKa/Ks");
            
            for (SyntenicBlock block : blocks) {
                for (SyntenicPair pair : block.getPairs()) {
                    String g1 = pair.getGeneId1();
                    String g2 = pair.getGeneId2();
                    
                    // Simple fuzzy match for .1 suffix
                    String s1 = seqs1.get(g1);
                    if (s1 == null) s1 = seqs1.get(g1 + ".1");
                    
                    String s2 = seqs2.get(g2);
                    if (s2 == null) s2 = seqs2.get(g2 + ".1");

                    if (s1 != null && s2 != null) {
                        // Align if lengths differ
                        String aligned1 = s1, aligned2 = s2;
                        if (s1.length() != s2.length()) {
                            String[] alignment = calculator.align(s1, s2);
                            aligned1 = alignment[0];
                            aligned2 = alignment[1];
                        }
                        
                        KaKsCalculator.Result res = calculator.calculate(aligned1, aligned2);
                        pw.printf(Locale.US, "%s\t%s\t%.6f\t%.6f\t%.6f%n", g1, g2, res.ka, res.ks, res.ratio);
                        count++;
                    }
                }
            }
        }

        System.out.println("[Phase 3/3] Done! Calculated Ka/Ks for " + count + " pairs.");
        System.out.println("Results saved to: " + output);
        
        return 0;
    }
}
