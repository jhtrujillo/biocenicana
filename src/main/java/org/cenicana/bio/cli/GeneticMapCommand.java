package org.cenicana.bio.cli;

import org.cenicana.bio.core.GeneticMapEngine;
import org.cenicana.bio.io.GeneticMapDashboardGenerator;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.concurrent.Callable;

@Command(name = "genetic-map",
         description = "Build genetic linkage maps (Linkage Groups & Marker Ordering) directly from VCF files.",
         mixinStandardHelpOptions = true,
         version = "genetic-map 1.0")
public class GeneticMapCommand implements Callable<Integer> {

    @Option(names = {"-i", "--input"}, description = "Path to the input VCF file containing population genotypes.", required = true)
    private String inputFile;

    @Option(names = {"-o", "--output"}, description = "Path to the output .map file to write genetic coordinates.", required = true)
    private String outputFile;

    @Option(names = {"--lod"}, description = "Minimum LOD score threshold to group markers into linkage groups.", defaultValue = "3.0")
    private double lod;

    @Option(names = {"--max-r"}, description = "Maximum recombination frequency (r) for linkage grouping.", defaultValue = "0.35")
    private double maxR;

    @Option(names = {"--mapping-function"}, description = "Mapping function to convert recombination to cM (kosambi, haldane).", defaultValue = "kosambi")
    private String mappingFunction;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level of the organism (e.g. 10 for sugarcane).", defaultValue = "10")
    private int ploidy;

    @Option(names = {"--single-dose"}, description = "Keep only single-dose (simplex x nulliplex) markers. Recommended for polyploid biparental populations.", defaultValue = "false")
    private boolean singleDose;

    @Option(names = {"--sd-chi2-p"}, description = "Minimum chi-square p-value for 1:1 segregation test (single-dose filter). Markers below this threshold are discarded.", defaultValue = "0.05")
    private double sdChi2P;

    @Option(names = {"--viz"}, description = "Path to generate an interactive HTML dashboard of the genetic map.")
    private String vizOutput;

    @Option(names = {"--thin-kb"}, description = "Physical thinning: keep 1 marker per window of N kb per chromosome. Reduces redundancy while preserving genome coverage. Recommended: 100-500 for sugarcane.", defaultValue = "0")
    private int thinKb;

    @Option(names = {"--pseudo-only"}, description = "Exclude contigs and scaffolds — keep only pseudochromosome markers.", defaultValue = "false")
    private boolean pseudoOnly;

    @Option(names = {"--min-lg-markers"}, description = "Minimum number of markers a Linkage Group must have to be included in the output. Default: 1 (all LGs).", defaultValue = "1")
    private int minLgMarkers;

    @Option(names = {"--genes-map"}, description = "TSV file with gene positions on the map (output of map_genes_to_map.py). Overlays candidate genes on the visual dashboard.")
    private String genesMapFile;

    @Override
    public Integer call() throws Exception {
        System.out.println("=================================================");
        System.out.println("BioJava: High-Performance Genetic Linkage Mapper");
        System.out.println("=================================================");
        System.out.println("Input VCF:          " + inputFile);
        System.out.println("Output Map:         " + outputFile);
        System.out.println("Min LOD Grouping:   " + lod);
        System.out.println("Max Recomb Limit:   " + maxR);
        System.out.println("Mapping Function:   " + mappingFunction.toUpperCase());
        System.out.println("Ploidy:             " + ploidy);
        System.out.println("Single-dose filter: " + (singleDose ? "YES (chi2 p >= " + sdChi2P + ")" : "NO"));
        System.out.println("Physical thinning:  " + (thinKb > 0 ? thinKb + " kb" : "NO"));
        System.out.println("Pseudo-chr only:    " + (pseudoOnly ? "YES" : "NO"));
        System.out.println("Min LG markers:     " + minLgMarkers);
        System.out.println("=================================================\n");

        long startTime = System.currentTimeMillis();

        GeneticMapEngine engine = new GeneticMapEngine(lod, maxR, mappingFunction);
        engine.setPloidy(ploidy);
        engine.setSingleDoseFilter(singleDose, sdChi2P);
        engine.setThinKb(thinKb);
        engine.setPseudoOnly(pseudoOnly);
        engine.setMinLgMarkers(minLgMarkers);
        engine.buildMap(inputFile, outputFile);

        if (vizOutput != null && !vizOutput.isEmpty()) {
            System.out.println("\n[MapViz] Generating interactive HTML dashboard...");
            GeneticMapDashboardGenerator.generate(outputFile, vizOutput, genesMapFile);
        }

        long endTime = System.currentTimeMillis();
        double elapsedSec = (endTime - startTime) / 1000.0;
        System.out.printf("\n🎉 Genetic map successfully built in %.2f seconds!\n", elapsedSec);

        return 0;
    }
}
