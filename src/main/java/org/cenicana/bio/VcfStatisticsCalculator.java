package org.cenicana.bio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import org.cenicana.bio.io.VcfFastReader;

public class VcfStatisticsCalculator {

	// Statistics
	public int numSnps = 0;
	public int numIndels = 0;
	public int numTransitions = 0;
	public int numTransversions = 0;

	// Per sample metrics
	public String[] sampleNames;
	public int[] sampleMissingCount;
	public long[] sampleTotalDepth;
	public int[] sampleDepthCount; // How many valid genotypes had depth info

	// Histograms
	// Site missingness bins: [0-0.1, 0.1-0.2, ... 0.9-1.0] (10 bins)
	public int[] siteMissingnessHistogram = new int[10];

	// Minor Allele Frequency (MAF) bins: [0-0.05, 0.05-0.10, ... 0.45-0.50] (10 bins)
	public int[] mafHistogram = new int[10];

	public void calculate(String vcfFile) throws IOException {
		sampleNames = VcfFastReader.getSampleIds(vcfFile);
		int numSamples = sampleNames.length;

		sampleMissingCount = new int[numSamples];
		sampleTotalDepth = new long[numSamples];
		sampleDepthCount = new int[numSamples];

		try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					continue;
				}

				String[] cols = line.split("\t");
				if (cols.length < 9 + numSamples) {
					continue; // Invalid row
				}

				String ref = cols[3];
				String alt = cols[4];
				
				// Identify SNP vs InDel
				if (ref.length() == 1 && alt.length() == 1 && !alt.equals(".")) {
					numSnps++;
					// Ts/Tv calculation
					if (isTransition(ref, alt)) {
						numTransitions++;
					} else {
						numTransversions++;
					}
				} else {
					numIndels++;
				}

				// Find DP/BSDP index in format field
				String format = cols[8];
				String[] formatFields = format.split(":");
				int depthIndex = -1;
				for (int i = 0; i < formatFields.length; i++) {
					if (formatFields[i].equals("BSDP") || formatFields[i].equals("DP")) {
						depthIndex = i;
						break;
					}
				}

				int siteMissing = 0;
				int totalAlleles = 0;
				int altAlleles = 0;

				// Parse genotypes
				for (int i = 0; i < numSamples; i++) {
					String genotypeData = cols[9 + i];
					
					if (genotypeData.startsWith("./.") || genotypeData.startsWith(".")) {
						sampleMissingCount[i]++;
						siteMissing++;
					} else {
						// Calculate MAF (very basic count based on 0/1 splits)
						String[] gtParts = genotypeData.split(":")[0].split("[/|]");
						for (String allele : gtParts) {
							if (!allele.equals(".")) {
								totalAlleles++;
								if (!allele.equals("0")) {
									altAlleles++;
								}
							}
						}

						// Calculate Depth
						if (depthIndex != -1) {
							String[] sampleFields = genotypeData.split(":");
							if (sampleFields.length > depthIndex && !sampleFields[depthIndex].equals(".")) {
								try {
									int depth = Integer.parseInt(sampleFields[depthIndex]);
									sampleTotalDepth[i] += depth;
									sampleDepthCount[i]++;
								} catch (NumberFormatException e) {
									// Ignore parsing error for depth
								}
							}
						}
					}
				}

				// Update Histograms
				double siteMissFreq = (double) siteMissing / numSamples;
				int binMissing = (int) (siteMissFreq * 10);
				if (binMissing >= 10) binMissing = 9;
				siteMissingnessHistogram[binMissing]++;

				if (totalAlleles > 0) {
					double altFreq = (double) altAlleles / totalAlleles;
					double maf = Math.min(altFreq, 1.0 - altFreq);
					int binMaf = (int) (maf * 20); // 0.0 to 0.5 -> 20 bins (size 0.025), wait 0 to 0.5 needs 10 bins of 0.05
					binMaf = (int) (maf / 0.05);
					if (binMaf >= 10) binMaf = 9;
					mafHistogram[binMaf]++;
				}
			}
		}
	}

	private boolean isTransition(String ref, String alt) {
		ref = ref.toUpperCase();
		alt = alt.toUpperCase();
		return (ref.equals("A") && alt.equals("G")) ||
				(ref.equals("G") && alt.equals("A")) ||
				(ref.equals("C") && alt.equals("T")) ||
				(ref.equals("T") && alt.equals("C"));
	}
}
