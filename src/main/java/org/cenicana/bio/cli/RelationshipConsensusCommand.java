package org.cenicana.bio.cli;

import org.cenicana.bio.core.PopulationStructureAnalyzer;
import org.cenicana.bio.core.PopulationStructureAnalyzer.PcaResult;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "rel-consensus",
         description = "Generates a consensus relationship report by cross-referencing PCA clusters, Genetic Distance, and Kinship.",
         mixinStandardHelpOptions = true,
         version = "rel-consensus 1.0")
public class RelationshipConsensusCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, description = "Path to the VCF file.", required = true)
    private String vcfFile;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level (e.g. 10 for sugarcane).", required = true)
    private int ploidy;

    @Option(names = {"-o", "--output"}, description = "Path to the output CSV file.", required = true)
    private String outputFile;

    @Option(names = {"--min-maf"}, description = "Minimum MAF.", defaultValue = "0.05")
    private double minMaf;

    @Option(names = {"--max-missing"}, description = "Maximum missingness.", defaultValue = "0.2")
    private double maxMissing;

    @Override
    public Integer call() throws Exception {
        System.out.println("=== BioJava: Relationship Consensus Report ===");
        
        PopulationStructureAnalyzer analyzer = new PopulationStructureAnalyzer();
        System.out.println("[Step 1/3] Running PCA, Distance and Kinship analysis...");
        PcaResult result = analyzer.computePCA(vcfFile, ploidy, 10, minMaf, maxMissing);

        System.out.println("[Step 2/3] Cross-referencing results and inferring relationships...");
        
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
            // Header
            pw.println("Sample1,Sample2,Same_PCA_Cluster,Distance_PC_Space,Kinship_VanRaden,Inferred_Relationship");

            int n = result.sampleNames.length;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    String s1 = result.sampleNames[i];
                    String s2 = result.sampleNames[j];
                    boolean sameCluster = result.clusterAssignments[i] == result.clusterAssignments[j];
                    double dist = result.distanceMatrix[i][j];
                    double kinship = result.kinshipMatrix[i][j];

                    String relationship = inferRelationship(kinship, dist);

                    pw.printf(Locale.US, "%s,%s,%s,%.4f,%.4f,%s\n", 
                        s1, s2, (sameCluster ? "YES" : "NO"), dist, kinship, relationship);
                }
            }
        }

        System.out.println("[Step 3/3] Report generated successfully: " + outputFile);
        return 0;
    }

    /**
     * Infers the biological relationship based on Kinship and Distance.
     * Thresholds are calibrated for polyploids but can be adjusted.
     */
    private String inferRelationship(double kinship, double dist) {
        if (kinship > 0.45) return "Duplicate/Clone";
        if (kinship > 0.20) return "Full-Sib/Parent-Offspring";
        if (kinship > 0.10) return "Half-Sib";
        if (kinship > 0.05) return "Distant Relative";
        if (kinship < -0.05) return "Unrelated (Divergent)";
        return "Unrelated";
    }
}
