package org.cenicana.bio.cli;

import java.io.File;
import java.util.concurrent.Callable;
import org.cenicana.bio.core.LinkageDisequilibriumCalculator;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command(name = "ld", description = "Calculate Pairwise Linkage Disequilibrium (r^2) using continuous Allele Dosages for Polyploids.")
public class LinkageDisequilibriumCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, required = true, description = "Path to the input VCF file.")
    private String vcfFile;

    @Option(names = {"-o", "--output"}, required = true, description = "Path to the output TSV file.")
    private String outputFile;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level of the organism. Used for semantic reference.", defaultValue = "2")
    private int ploidy;

    @Option(names = {"-m", "--min-maf"}, description = "Minimum MAF. Markers with lower MAF are ignored to prevent artificial r^2 inflation.", defaultValue = "0.05")
    private double minMaf;

    @Option(names = {"-w", "--window-size"}, description = "Maximum distance (in bp) between two markers to calculate LD.", defaultValue = "50000")
    private int windowSize;

    @Option(names = {"--html"}, description = "Generate an interactive HTML dashboard plotting the LD Decay curve.")
    private boolean generateHtml;

    @Option(names = {"--bin-size"}, description = "Bin size in bp for the LD decay curve graph.", defaultValue = "1000")
    private int binSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads for parallel processing.", defaultValue = "-1")
    private int threads;

    @Override
    public Integer call() throws Exception {
        File f = new File(vcfFile);
        if (!f.exists()) {
            System.err.println("Error: VCF file not found: " + vcfFile);
            return 1;
        }

        System.out.println("=================================================");
        System.out.println("BioJava Polyploid LD Calculator");
        System.out.println("=================================================");
        System.out.println("Input:  " + vcfFile);
        System.out.println("Output: " + outputFile);
        System.out.println("Parameters:");
        System.out.println(" - Ploidy assumption: " + ploidy);
        System.out.println(" - Min MAF: " + minMaf);
        System.out.println(" - Sliding Window (bp): " + windowSize);
        if (generateHtml) System.out.println(" - HTML Dashboard: yes (bin size: " + binSize + " bp)");
        System.out.println("=================================================\n");

        LinkageDisequilibriumCalculator ldCalc = new LinkageDisequilibriumCalculator();
        ldCalc.setPloidy(ploidy);
        ldCalc.setMinMaf(minMaf);
        ldCalc.setMaxDistanceBp(windowSize);
        ldCalc.setGenerateHtml(generateHtml);
        ldCalc.setBinSizeBp(binSize);
        
        int numThreads = (threads <= 0) ? Runtime.getRuntime().availableProcessors() : threads;
        ldCalc.computeLD(vcfFile, outputFile, numThreads);

        return 0;
    }
}
