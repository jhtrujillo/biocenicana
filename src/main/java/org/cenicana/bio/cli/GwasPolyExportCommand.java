package org.cenicana.bio.cli;

import org.cenicana.bio.core.GwasPolyExporter;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.io.File;
import java.util.concurrent.Callable;

@Command(name = "gwaspoly-export", description = "Export VCF to GWASpoly ACGT format (ploidy-aware strings)")
public class GwasPolyExportCommand implements Callable<Integer> {

    @Option(names = {"-v", "--vcf"}, required = true, description = "Input VCF file")
    private String vcfFile;

    @Option(names = {"-o", "--output"}, required = true, description = "Output CSV file (GWASpoly format)")
    private String outputFile;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level (e.g., 10 for sugarcane)", defaultValue = "10")
    private int ploidy;

    @Override
    public Integer call() throws Exception {
        File f = new File(vcfFile);
        if (!f.exists()) {
            System.err.println("Error: VCF file not found: " + vcfFile);
            return 1;
        }

        System.out.println("=================================================");
        System.out.println("BioCenicana GWASpoly Exporter");
        System.out.println("=================================================");
        System.out.println("Input:  " + vcfFile);
        System.out.println("Output: " + outputFile);
        System.out.println("Ploidy: " + ploidy);
        System.out.println("=================================================\n");

        GwasPolyExporter exporter = new GwasPolyExporter();
        exporter.exportToACGT(vcfFile, outputFile, ploidy);

        System.out.println("\n✅  Export successful.");
        return 0;
    }
}
