package org.cenicana.bio.cli;

import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.PhylogenyTreeBuilder;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "snp-tree",
         description = "Builds a SNP-based phylogeny (Neighbor-Joining tree) from a VCF file.",
         mixinStandardHelpOptions = true,
         version = "snp-tree 1.0")
public class SnpTreeCommand implements Callable<Integer> {

    enum DistanceMethod {
        MANHATTAN("manhattan"),
        EUCLIDEAN("euclidean"),
        NEI("nei"),
        ROGERS("rogers"),
        IBS("ibs"),
        P_DISTANCE("p-distance");

        private final String stringValue;
        DistanceMethod(String stringValue) { this.stringValue = stringValue; }
        @Override public String toString() { return stringValue; }
    }

    static class MethodConverter implements picocli.CommandLine.ITypeConverter<DistanceMethod> {
        @Override
        public DistanceMethod convert(String value) throws Exception {
            for (DistanceMethod m : DistanceMethod.values()) {
                if (m.stringValue.equalsIgnoreCase(value)) return m;
            }
            throw new Exception("Invalid method: " + value);
        }
    }

    @Option(names = {"-v", "--vcf"}, description = "Path to the VCF file.", required = true)
    private String vcfFile;

    @Option(names = {"-o", "--output"}, description = "Path to the output Newick file.", required = true)
    private String outputFile;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level (e.g. 10 for sugarcane).", required = true)
    private int ploidy;

    @Option(names = {"-c", "--caller"}, description = "Variant caller (auto, ngsep, gatk, freebayes).", defaultValue = "auto")
    private String caller;

    @Option(names = {"-md", "--min-depth"}, description = "Minimum depth for genotype trust.", defaultValue = "0")
    private int minDepth;

    @Option(names = {"-t", "--threads"}, description = "Number of threads.", defaultValue = "-1")
    private int threads;

    @Option(names = {"--pca"}, description = "Optional: path to a PCA CSV file (from pop-structure) to color nodes by cluster.")
    private String pcaFile;

    @Option(names = {"-m", "--method"},
            description = "Genetic distance method (manhattan, euclidean, nei, rogers, ibs, p-distance).",
            defaultValue = "manhattan",
            converter = MethodConverter.class)
    private DistanceMethod method;

    @Override
    public Integer call() throws Exception {
        System.out.println("=== BioJava: SNP-based Phylogeny (Neighbor-Joining) ===");

        int numThreads = (threads <= 0) ? Runtime.getRuntime().availableProcessors() : threads;
        AlleleDosageCalculator calculator = new AlleleDosageCalculator();

        System.out.println("[Step 1/2] Computing genetic distance matrix from VCF using " + method + " method...");
        AlleleDosageCalculator.DistanceResult dr = calculator.computeDistanceMatrix(vcfFile, ploidy, caller, minDepth, false, numThreads, method.toString());

        System.out.println("[Step 2/2] Building Neighbor-Joining tree...");
        PhylogenyTreeBuilder treeBuilder = new PhylogenyTreeBuilder();
        String newick = treeBuilder.buildNewick(dr.matrix, dr.sampleIds);

        // Save Newick
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
            pw.println(newick);
        }

        // Parse PCA clusters if provided
        String pcaJson = "null";
        if (pcaFile != null && !pcaFile.isBlank()) {
            pcaJson = parsePcaClusters(pcaFile);
            System.out.println("[Info] PCA cluster data loaded from: " + pcaFile);
        }

        // Save HTML Report
        String htmlFile = outputFile.endsWith(".nwk") ? outputFile.replace(".nwk", ".html") : outputFile + ".html";
        generateHtmlReport(newick, pcaJson, htmlFile);

        System.out.println("[Done] Tree saved in Newick format to: " + outputFile);
        System.out.println("[Done] Interactive viewer saved to: " + htmlFile);
        System.out.println("You can open the HTML file in any browser to explore the tree.");

        return 0;
    }

    /**
     * Parses the PCA CSV and returns a JSON object mapping sampleId -> {kmeans, dbscan, gmm}
     */
    private String parsePcaClusters(String csvPath) throws Exception {
        StringBuilder sb = new StringBuilder("{");
        boolean first = true;
        List<String> lines = Files.readAllLines(Paths.get(csvPath));
        if (lines.isEmpty()) return "null";

        // Parse header to find column indices
        String[] header = lines.get(0).split(",");
        int idxSample = -1, idxKMeans = -1, idxDBSCAN = -1, idxGMM = -1;
        for (int i = 0; i < header.length; i++) {
            String h = header[i].trim();
            if (h.equalsIgnoreCase("Sample"))          idxSample = i;
            else if (h.equalsIgnoreCase("KMeans_Cluster")) idxKMeans = i;
            else if (h.equalsIgnoreCase("DBSCAN_Cluster")) idxDBSCAN = i;
            else if (h.equalsIgnoreCase("GMM_Cluster"))    idxGMM    = i;
        }
        if (idxSample < 0) return "null";

        for (int r = 1; r < lines.size(); r++) {
            String line = lines.get(r).trim();
            if (line.isEmpty()) continue;
            String[] cols = line.split(",", -1);
            String sample = cols[idxSample].trim();
            String kmeans = (idxKMeans >= 0 && idxKMeans < cols.length) ? cols[idxKMeans].trim() : "?";
            String dbscan = (idxDBSCAN >= 0 && idxDBSCAN < cols.length) ? cols[idxDBSCAN].trim() : "?";
            String gmm    = (idxGMM    >= 0 && idxGMM    < cols.length) ? cols[idxGMM].trim()    : "?";

            if (!first) sb.append(",");
            sb.append("\"").append(sample).append("\":{")
              .append("\"kmeans\":\"").append(kmeans).append("\",")
              .append("\"dbscan\":\"").append(dbscan).append("\",")
              .append("\"gmm\":\"").append(gmm).append("\"}");
            first = false;
        }
        sb.append("}");
        return sb.toString();
    }

    private void generateHtmlReport(String newick, String pcaJson, String outputFile) throws Exception {
        InputStream is = getClass().getResourceAsStream("/phylogeny_template.html");
        if (is == null) {
            System.err.println("[Warning] Could not find phylogeny_template.html in resources.");
            return;
        }
        Scanner s = new Scanner(is).useDelimiter("\\A");
        String template = s.hasNext() ? s.next() : "";

        String html = template
                .replace("/*NEWICK_DATA*/", newick)
                .replace("/*PCA_DATA*/", pcaJson);

        try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
            pw.println(html);
        }
    }
}
