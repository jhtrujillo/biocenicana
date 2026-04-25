package org.cenicana.bio.core;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.cenicana.bio.io.VcfFastReader;
import org.cenicana.bio.utils.HweUtils;

/**
 * A fast, streaming VCF filter that mimics NGSEP's VCFFilter.
 * Can filter by MAF, missingness, HWE p-value, and biallelic SNPs only.
 */
public class VcfFilter {

	private double minMaf = 0.0;
	private double maxMissingness = 1.0;
	private double minHwePValue = 0.0;
	private boolean onlyBiallelicSnps = false;
	private int ploidy = 2;
	private double minEh = 0.0;
	private int topN = -1;
	private int minDepth = 5;
	private boolean strictGenotypes = false; // If true, discard ./. even if counts exist

	// Helper class for PriorityQueue to keep Top N polymorphic markers
	private static class VcfLineScore implements Comparable<VcfLineScore> {
		long order;
		String line;
		double eh;
		public VcfLineScore(long order, String line, double eh) {
			this.order = order;
			this.line = line;
			this.eh = eh;
		}
		@Override
		public int compareTo(VcfLineScore o) {
			return Double.compare(this.eh, o.eh); // Min-heap (lowest EH at the head)
		}
	}

	public void setMinMaf(double minMaf) { this.minMaf = minMaf; }
	public void setMaxMissingness(double maxMissingness) { this.maxMissingness = maxMissingness; }
	public void setMinHwePValue(double minHwePValue) { this.minHwePValue = minHwePValue; }
	public void setOnlyBiallelicSnps(boolean onlyBiallelicSnps) { this.onlyBiallelicSnps = onlyBiallelicSnps; }
	public void setPloidy(int ploidy) { this.ploidy = ploidy; }
	public void setMinEh(double minEh) { this.minEh = minEh; }
	public void setTopN(int topN) { this.topN = topN; }
	public void setMinDepth(int minDepth) { this.minDepth = minDepth; }
	public void setStrictGenotypes(boolean strictGenotypes) { this.strictGenotypes = strictGenotypes; }

	public void filter(String inputVcf, String outputVcf, int threads) throws IOException {
		String[] sampleNames = VcfFastReader.getSampleIds(inputVcf);
		int numSamples = sampleNames.length;

		AtomicInteger keptVariants = new AtomicInteger(0);
		AtomicInteger totalVariants = new AtomicInteger(0);
		
		// If Top N is active, we use the original sequential implementation 
		// because Top N requires global knowledge and sorting by EH.
		if (topN > 0) {
			filterSequential(inputVcf, outputVcf, numSamples);
			return;
		}

		ExecutorService executor = Executors.newFixedThreadPool(threads);
		final Map<Long, String> resultBuffer = Collections.synchronizedMap(new TreeMap<>());
		final AtomicLong nextToWrite = new AtomicLong(1);
		final Object lock = new Object();

		try (BufferedReader br = new BufferedReader(new FileReader(inputVcf));
			 PrintWriter pw = new PrintWriter(new FileWriter(outputVcf))) {

			String line;
			long sequenceCounter = 0;
			
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					if (line.startsWith("#CHROM")) {
						pw.println("##BioCenicana_VCFFilter=\"minMaf=" + minMaf + ",maxMissingness=" + maxMissingness + ",minHwePValue=" + minHwePValue + ",onlyBiallelicSnps=" + onlyBiallelicSnps + ",ploidy=" + ploidy + ",minEh=" + minEh + ",topN=" + topN + ",threads=" + threads + "\"");
					}
					pw.println(line);
					continue;
				}

				totalVariants.incrementAndGet();
				sequenceCounter++;
				final long currentSeq = sequenceCounter;
				final String currentLine = line;

				executor.submit(() -> {
					String result = processLine(currentLine, numSamples);
					
					synchronized (lock) {
						resultBuffer.put(currentSeq, result);
						while (resultBuffer.containsKey(nextToWrite.get())) {
							String toWrite = resultBuffer.remove(nextToWrite.get());
							if (toWrite != null) {
								pw.println(toWrite);
								keptVariants.incrementAndGet();
							}
							nextToWrite.incrementAndGet();
						}
						lock.notifyAll();
					}
				});

				// Memory protection: if workers are too fast reading but slow writing
				if (sequenceCounter - nextToWrite.get() > 10000) {
					synchronized (lock) {
						while (sequenceCounter - nextToWrite.get() > 2000) {
							try { lock.wait(50); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
						}
					}
				}
			}

			executor.shutdown();
			try {
				executor.awaitTermination(1, TimeUnit.HOURS);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}
		}

		System.out.println("[VcfFilter] Processing complete (Parallel).");
		System.out.println("[VcfFilter] Total variants processed: " + totalVariants.get());
		System.out.println("[VcfFilter] Variants kept: " + keptVariants.get());
		System.out.println("[VcfFilter] Variants removed: " + (totalVariants.get() - keptVariants.get()));
	}

	private void filterSequential(String inputVcf, String outputVcf, int numSamples) throws IOException {
		int totalVariants = 0;
		int keptVariants = 0;
		PriorityQueue<VcfLineScore> topQueue = new PriorityQueue<>();

		try (BufferedReader br = new BufferedReader(new FileReader(inputVcf));
			 PrintWriter pw = new PrintWriter(new FileWriter(outputVcf))) {

			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					if (line.startsWith("#CHROM")) {
						pw.println("##BioCenicana_VCFFilter=\"minMaf=" + minMaf + ",maxMissingness=" + maxMissingness + ",minHwePValue=" + minHwePValue + ",onlyBiallelicSnps=" + onlyBiallelicSnps + ",ploidy=" + ploidy + ",minEh=" + minEh + ",topN=" + topN + ",mode=sequential\"");
					}
					pw.println(line);
					continue;
				}

				totalVariants++;
				String result = processLine(line, numSamples);
				if (result != null) {
					if (topN > 0) {
						double eh = calculateEh(line, numSamples);
						topQueue.add(new VcfLineScore(totalVariants, line, eh));
						if (topQueue.size() > topN) topQueue.poll();
					} else {
						pw.println(line);
						keptVariants++;
					}
				}
			}

			if (topN > 0) {
				List<VcfLineScore> sortedList = new ArrayList<>(topQueue);
				sortedList.sort(Comparator.comparingLong(s -> s.order));
				for (VcfLineScore s : sortedList) {
					pw.println(s.line);
					keptVariants++;
				}
			}
		}
		System.out.println("[VcfFilter] Processing complete (Sequential/TopN).");
		System.out.println("[VcfFilter] Total variants processed: " + totalVariants);
		System.out.println("[VcfFilter] Variants kept: " + keptVariants);
	}

	private String processLine(String line, int numSamples) {
		String[] cols = line.split("\t");
		if (cols.length < 9 + numSamples) return null;

		String ref = cols[3];
		String alt = cols[4];
		String[] altAllelesList = alt.equals(".") ? new String[0] : alt.split(",");
		int numAllelesAtSite = 1 + altAllelesList.length;

		boolean isBiallelicSnp = ref.length() == 1
			&& altAllelesList.length == 1
			&& altAllelesList[0].length() == 1
			&& !alt.equals(".");

		if (onlyBiallelicSnps && !isBiallelicSnp) return null;

		int siteMissing = 0;
		int totalGenotypedAlleles = 0;
		int[] siteAlleleCounts = new int[numAllelesAtSite];
		int n00 = 0, n01 = 0, n11 = 0; 
		int[] polyDosageCounts = (ploidy > 2 && isBiallelicSnp) ? new int[ploidy + 1] : null;
		int nGenotypedPoly = 0;

		int fmtAdIdx = -1, fmtBsdpIdx = -1, fmtRoIdx = -1, fmtAoIdx = -1;
		String[] formatFields = cols[8].split(":");
		for (int f = 0; f < formatFields.length; f++) {
			switch (formatFields[f]) {
				case "AD":   fmtAdIdx   = f; break;
				case "BSDP": fmtBsdpIdx = f; break;
				case "RO":   fmtRoIdx   = f; break;
				case "AO":   fmtAoIdx   = f; break;
			}
		}

		for (int i = 0; i < numSamples; i++) {
			String[] gData = cols[9 + i].split(":");
			String gt = gData[0];

			if (ploidy > 2 && isBiallelicSnp) {
				if (strictGenotypes && gt.startsWith(".")) {
					siteMissing++;
					continue;
				}
				double[] counts = AlleleDosageCalculator.getRefAltCounts(gData, fmtAdIdx, fmtBsdpIdx, fmtRoIdx, fmtAoIdx, ref, alt);
				if (counts != null && counts[0] + counts[1] >= minDepth) {
					double refReads = counts[0];
					double altReads = counts[1];
					double totalReads = refReads + altReads;
					int dosage = (int) Math.round((altReads / totalReads) * ploidy);
					if (dosage >= 0 && dosage <= ploidy) {
						polyDosageCounts[dosage]++;
						siteAlleleCounts[0] += (ploidy - dosage);
						siteAlleleCounts[1] += dosage;
						totalGenotypedAlleles += ploidy;
						nGenotypedPoly++;
					} else {
						siteMissing++;
					}
				} else {
					siteMissing++;
				}
			} else {
				if (gt.equals(".") || gt.startsWith("./.") || gt.startsWith("./")) {
					siteMissing++;
				} else {
					String[] alleles = gt.split("[/|]");
					boolean valid = true;
					for (String a : alleles) {
						if (a.equals(".")) { valid = false; break; }
						int idx = Integer.parseInt(a);
						if (idx < numAllelesAtSite) {
							siteAlleleCounts[idx]++;
							totalGenotypedAlleles++;
						} else { valid = false; }
					}
					if (valid && isBiallelicSnp && alleles.length == 2) {
						if (alleles[0].equals("0") && alleles[1].equals("0")) n00++;
						else if (alleles[0].equals("1") && alleles[1].equals("1")) n11++;
						else n01++;
					}
					if (!valid) siteMissing++;
				}
			}
		}

		if ((double) siteMissing / numSamples > maxMissingness) return null;

		if (totalGenotypedAlleles > 0) {
			int minAc = Integer.MAX_VALUE;
			for (int c : siteAlleleCounts) if (c > 0 && c < minAc) minAc = c;
			if ((double) minAc / totalGenotypedAlleles < minMaf) return null;

			if (isBiallelicSnp && minHwePValue > 0.0) {
				if (ploidy > 2 && polyDosageCounts != null && nGenotypedPoly >= 5) {
					double totalAltDosage = 0;
					for (int d = 0; d <= ploidy; d++) totalAltDosage += polyDosageCounts[d] * d;
					double q = totalAltDosage / (nGenotypedPoly * ploidy);
					double p = 1.0 - q;
					double chi2 = 0.0;
					double[] logFacs = new double[ploidy + 1];
					for(int j=1; j<=ploidy; j++) logFacs[j] = logFacs[j-1] + Math.log(j);
					for (int d = 0; d <= ploidy; d++) {
						double logComb = logFacs[ploidy] - logFacs[d] - logFacs[ploidy - d];
						double eProb = Math.exp(logComb + (ploidy - d)*Math.log(Math.max(1e-10, p)) + d*Math.log(Math.max(1e-10, q)));
						double eCount = Math.max(1e-5, nGenotypedPoly * eProb);
						chi2 += Math.pow(polyDosageCounts[d] - eCount, 2) / eCount;
					}
					if (chiSquarePValue(chi2, ploidy - 1) < minHwePValue) return null;
				} else if (ploidy <= 2 && (n00+n01+n11) >= 5) {
					if (org.cenicana.bio.utils.HweUtils.calculateHweFisherExactTest(n00, n01, n11) < minHwePValue) return null;
				}
			}
			double eh = 0.0;
			for (int c : siteAlleleCounts) if (c > 0) eh += Math.pow((double) c / totalGenotypedAlleles, 2);
			if (1.0 - eh < minEh) return null;
		} else return null;

		return line;
	}

	private double calculateEh(String line, int numSamples) {
		// Helper for Top N mode
		String[] cols = line.split("\t");
		String ref = cols[3];
		String alt = cols[4];
		String[] altAllelesList = alt.equals(".") ? new String[0] : alt.split(",");
		int numAllelesAtSite = 1 + altAllelesList.length;
		int totalGenotypedAlleles = 0;
		int[] siteAlleleCounts = new int[numAllelesAtSite];

		for (int i = 0; i < numSamples; i++) {
			String gt = cols[9 + i].split(":")[0];
			if (!gt.startsWith(".")) {
				for (String a : gt.split("[/|]")) {
					if (!a.equals(".")) {
						int idx = Integer.parseInt(a);
						if (idx < numAllelesAtSite) { siteAlleleCounts[idx]++; totalGenotypedAlleles++; }
					}
				}
			}
		}
		double eh = 0.0;
		if (totalGenotypedAlleles > 0) {
			for (int c : siteAlleleCounts) eh += Math.pow((double) c / totalGenotypedAlleles, 2);
		}
		return 1.0 - eh;
	}


	private double chiSquarePValue(double x, int df) {
		if (x <= 0 || df < 1) return 1.0;
		double dfD = df;
		// Wilson-Hilferty transformation for Chi-Square to Normal distribution
		double Z = (Math.pow(x / dfD, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * dfD))) / Math.sqrt(2.0 / (9.0 * dfD));
		return 0.5 * erfc(Z / Math.sqrt(2.0));
	}

	private double erfc(double x) {
		double z = Math.abs(x);
		double t = 1.0 / (1.0 + 0.5 * z);
		double ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
		return x >= 0 ? ans : 2.0 - ans;
	}
}
