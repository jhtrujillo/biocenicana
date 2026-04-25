package org.cenicana.bio.cli;

import org.cenicana.bio.core.ComparativeGenomicsAnalyzer;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.util.concurrent.Callable;

/**
 * Command for Comparative Genomics analysis.
 */
@Command(name = "comp-gen", 
         mixinStandardHelpOptions = true, 
         description = "Integrate GFF, CDS, Proteins and Synteny results (McScanX/Synmap2).")
public class ComparativeGenomicsCommand implements Callable<Integer> {

    @Option(names = {"--gff1"}, required = true, description = "GFF file for Genome 1.")
    private String gff1;

    @Option(names = {"--gff2"}, required = true, description = "GFF file for Genome 2.")
    private String gff2;

    @Option(names = {"--collinearity"}, required = true, description = "McScanX .collinearity file.")
    private String collinearity;

    @Option(names = {"--cds1"}, description = "CDS FASTA file for Genome 1.")
    private String cds1;

    @Option(names = {"--cds2"}, description = "CDS FASTA file for Genome 2.")
    private String cds2;

    @Option(names = {"--prot1"}, description = "Protein FASTA file for Genome 1.")
    private String prot1;

    @Option(names = {"--prot2"}, description = "Protein FASTA file for Genome 2.")
    private String prot2;

    @Option(names = {"-o", "--output"}, defaultValue = "comparative_report.tsv", description = "Output TSV file.")
    private String output;

    @Override
    public Integer call() throws Exception {
        ComparativeGenomicsAnalyzer analyzer = new ComparativeGenomicsAnalyzer();
        analyzer.runAnalysis(gff1, gff2, collinearity, cds1, cds2, prot1, prot2, output);
        return 0;
    }
}
