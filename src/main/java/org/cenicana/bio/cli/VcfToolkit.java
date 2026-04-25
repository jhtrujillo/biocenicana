package org.cenicana.bio.cli;

import picocli.CommandLine;
import picocli.CommandLine.Command;

/**
 * Main entrance for the BioCenicana Bioinformatics Toolkit.
 * This class aggregates all specialized subcommands for VCF processing and analysis.
 */
@Command(name = "biocenicana", 
         mixinStandardHelpOptions = true, 
         version = "1.0", 
         description = "Bioinformatics tools for Cenicana (Genomics of Polyploids)", 
         subcommands = {
             AlleleDosageCommand.class,
             GeneticDistanceCommand.class,
             VcfStatsCommand.class,
             VcfFilterCommand.class,
             LinkageDisequilibriumCommand.class,
             GwasPolyExportCommand.class,
             PopStructureCommand.class,
             JoinMapCommand.class,
             SnpExplorerCommand.class,
             VcfMergeCommand.class,
             ComparativeGenomicsCommand.class
         })
public class VcfToolkit {

    public static void main(String[] args) {
        int exitCode = new CommandLine(new VcfToolkit()).execute(args);
        if (exitCode != 0) {
            System.exit(exitCode);
        }
    }
}
