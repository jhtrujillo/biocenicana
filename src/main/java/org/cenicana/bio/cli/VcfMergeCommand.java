package org.cenicana.bio.cli;

import org.cenicana.bio.core.VcfMerger;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;

@Command(name = "vcf-merge", 
         mixinStandardHelpOptions = true, 
         description = "Merges multiple VCF files into one, handling union of samples and filling missing data.")
public class VcfMergeCommand implements Callable<Integer> {

    @Option(names = { "-i", "--input" }, required = true, split = ",", 
            description = "Comma-separated list of VCF files to merge.")
    private String[] inputFiles;

    @Option(names = { "-o", "--output" }, required = true, 
            description = "Path to the output merged VCF file.")
    private String outputFile;

    @Override
    public Integer call() throws Exception {
        System.out.println("🧬 BioJava - VCF Merge");
        System.out.println("-------------------------------------");

        if (inputFiles.length < 2) {
            System.err.println("Error: At least two input files are required for merging.");
            return 1;
        }

        List<String> files = Arrays.asList(inputFiles);
        VcfMerger merger = new VcfMerger(files, outputFile);

        long startTime = System.currentTimeMillis();
        try {
            merger.merge();
            long endTime = System.currentTimeMillis();
            System.out.println("\n✅  Merge completed successfully!");
            System.out.println("Output: " + outputFile);
            System.out.println("Time: " + (endTime - startTime) / 1000.0 + "s");
        } catch (Exception e) {
            System.err.println("\n❌ Error during merge: " + e.getMessage());
            e.printStackTrace();
            return 1;
        }

        return 0;
    }
}
