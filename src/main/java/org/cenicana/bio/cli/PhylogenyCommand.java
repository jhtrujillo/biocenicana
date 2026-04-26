package org.cenicana.bio.cli;

import org.cenicana.bio.core.PhylogenyAnalyzer;
import org.cenicana.bio.io.CollinearityParser;
import org.cenicana.bio.io.GffParser;
import org.cenicana.bio.model.Gene;
import org.cenicana.bio.model.SyntenicBlock;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

@Command(name = "phylo", mixinStandardHelpOptions = true, version = "1.0",
        description = "Phylogeny Module: Calculates evolutionary distances between orthologs.")
public class PhylogenyCommand implements Callable<Integer> {

    @Option(names = {"--collinearity"}, description = "Collinearity file from McScanX", required = true)
    private String collinearity;

    @Option(names = {"--gff1"}, description = "GFF file for Genome 1", required = true)
    private String gff1;

    @Option(names = {"--gff2"}, description = "GFF file for Genome 2", required = true)
    private String gff2;

    @Option(names = {"--fasta1"}, description = "FASTA file (CDS/Protein) for Genome 1", required = true)
    private String fasta1;

    @Option(names = {"--fasta2"}, description = "FASTA file (CDS/Protein) for Genome 2", required = true)
    private String fasta2;

    @Option(names = {"-o", "--output"}, description = "Output prefix", defaultValue = "phylogeny_results")
    private String output;

    @Override
    public Integer call() throws Exception {
        System.out.println("=== BioCenicana: Phylogeny Module ===");
        
        GffParser gffParser = new GffParser();
        Map<String, Gene> map1 = loadGenes(gffParser.parse(gff1));
        Map<String, Gene> map2 = loadGenes(gffParser.parse(gff2));

        CollinearityParser colParser = new CollinearityParser();
        List<SyntenicBlock> blocks = colParser.parse(collinearity);

        PhylogenyAnalyzer analyzer = new PhylogenyAnalyzer();
        analyzer.runPhylogeny(blocks, map1, map2, fasta1, fasta2, output);

        System.out.println("=== Phylogeny Analysis Complete ===");
        return 0;
    }

    private Map<String, Gene> loadGenes(List<Gene> list) {
        Map<String, Gene> map = new HashMap<>();
        for (Gene g : list) map.put(g.getId(), g);
        return map;
    }
}
