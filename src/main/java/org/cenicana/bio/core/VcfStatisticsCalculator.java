package org.cenicana.bio.core;
import org.cenicana.bio.utils.FileUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.cenicana.bio.io.VcfFastReader;
import org.cenicana.bio.utils.HweUtils;

/**
 * Streaming VCF statistics calculator with multi-threading support.
 */
public class VcfStatisticsCalculator {

	// ── Global variant counts ────────────────────────────────────────────────
	public int numSnps          = 0;
	public int numIndels        = 0;
	public int numTransitions   = 0;
	public int numTransversions = 0;

	// ── Per-sample metrics ───────────────────────────────────────────────────
	public String[] sampleNames;
	public int[]    sampleMissingCount;
	public int[]    sampleGenotypedCount;
	public int[]    sampleHetCount;
	public int[]    sampleHomoRefCount;
	public int[]    sampleHomoAltCount;
	public int[]    sampleNonRefCount;
	public int[]    sampleTsCount;
	public int[]    sampleTvCount;
	public int[]    sampleRareAlleleCount;
	public long[]   sampleTotalDepth;
	public int[]    sampleDepthCount;

	public double[] sampleObsHet;
	public double[] sampleFstat;

	private double ehSum   = 0;
	private int    ehCount = 0;
	public  double meanEH  = 0;

	public Map<String, Integer> variantsPerChromosome = new LinkedHashMap<>();
	public static final int DENSITY_BIN_SIZE = 1000000;
	public Map<String, Map<Integer, Integer>> binnedDensity = new LinkedHashMap<>();

	public int[] siteMissingnessHistogram = new int[10];
	public int[] mafHistogram             = new int[50];
	public int[] ehHistogram              = new int[10];

	public int[]  hweChiSqHistogram  = new int[5];
	public int    numHweViolations   = 0;
	public int    numHweTested       = 0;
	public double meanChiSq          = 0;
	private double chiSqSum          = 0;

	public int numMonomorphic  = 0;
	public int numBiallelic    = 0;
	public int numMultiallelic = 0;

	public double tajimaD       = Double.NaN;
	public int    numSegSites   = 0;
	public double piHat         = 0;
	public double thetaW        = 0;
	private double piPerSiteSum = 0;
	private long   totalAllelesForTajima = 0;
	private int    tajimaSiteCount       = 0;

	private Map<String, Integer> samplePopIndex = null;
	public  String[]             populationNames = null;
	public  double[][]           pairwiseFst    = null;

	public void loadPopulationMap(String popFile) throws IOException {
		Map<String, String> sampleToPop = new HashMap<>();
		List<String> popOrder = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(popFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) continue;
				String[] parts = line.split("\t");
				if (parts.length < 2) continue;
				String sample = parts[0].trim();
				String pop    = parts[1].trim();
				sampleToPop.put(sample, pop);
				if (!popOrder.contains(pop)) popOrder.add(pop);
			}
		}
		populationNames  = popOrder.toArray(new String[0]);
		samplePopIndex   = new HashMap<>();
		for (int i = 0; i < populationNames.length; i++) {
			final int idx = i;
			sampleToPop.forEach((s, p) -> {
				if (p.equals(populationNames[idx])) samplePopIndex.put(s, idx);
			});
		}
	}

	public void calculate(String vcfFile, int threads) throws IOException {
		sampleNames = VcfFastReader.getSampleIds(vcfFile);
		int numSamples = sampleNames.length;

		sampleMissingCount    = new int[numSamples];
		sampleGenotypedCount  = new int[numSamples];
		sampleHetCount        = new int[numSamples];
		sampleHomoRefCount    = new int[numSamples];
		sampleHomoAltCount    = new int[numSamples];
		sampleNonRefCount     = new int[numSamples];
		sampleTsCount         = new int[numSamples];
		sampleTvCount         = new int[numSamples];
		sampleRareAlleleCount = new int[numSamples];
		sampleTotalDepth      = new long[numSamples];
		sampleDepthCount      = new int[numSamples];

		final int[] samplePopIdxArr;
		final int   numPops;
		if (samplePopIndex != null && populationNames != null) {
			numPops      = populationNames.length;
			samplePopIdxArr = new int[numSamples];
			for (int i = 0; i < numSamples; i++) {
				Integer pidx = samplePopIndex.get(sampleNames[i]);
				samplePopIdxArr[i] = (pidx != null) ? pidx : -1;
			}
		} else {
			numPops = 0;
			samplePopIdxArr = null;
		}

		ExecutorService executor = Executors.newFixedThreadPool(threads);
		List<LocalStats> localStatsList = Collections.synchronizedList(new ArrayList<>());

		try (BufferedReader br = new BufferedReader(new FileReader(vcfFile))) {
			String line;
			List<String> chunk = new ArrayList<>(1000);
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) continue;
				chunk.add(line);
				if (chunk.size() >= 1000) {
					final List<String> currentChunk = chunk;
					executor.submit(() -> {
						LocalStats ls = new LocalStats(numSamples, numPops);
						for (String l : currentChunk) {
							processLine(l, numSamples, ls, samplePopIdxArr);
						}
						localStatsList.add(ls);
					});
					chunk = new ArrayList<>(1000);
				}
			}
			if (!chunk.isEmpty()) {
				final List<String> lastChunk = chunk;
				executor.submit(() -> {
					LocalStats ls = new LocalStats(numSamples, numPops);
					for (String l : lastChunk) {
						processLine(l, numSamples, ls, samplePopIdxArr);
					}
					localStatsList.add(ls);
				});
			}
			executor.shutdown();
			executor.awaitTermination(1, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		}

		double[][][] fstAccum = (numPops > 0) ? new double[numPops][numPops][3] : null;

		for (LocalStats ls : localStatsList) {
			numSnps += ls.numSnps; numIndels += ls.numIndels; numTransitions += ls.numTransitions; numTransversions += ls.numTransversions;
			ehSum += ls.ehSum; ehCount += ls.ehCount; numHweViolations += ls.numHweViolations; numHweTested += ls.numHweTested;
			chiSqSum += ls.chiSqSum; numMonomorphic += ls.numMonomorphic; numBiallelic += ls.numBiallelic; numMultiallelic += ls.numMultiallelic;
			numSegSites += ls.numSegSites; piPerSiteSum += ls.piPerSiteSum; totalAllelesForTajima += ls.totalAllelesForTajima; tajimaSiteCount += ls.tajimaSiteCount;

			for (int i = 0; i < numSamples; i++) {
				sampleMissingCount[i] += ls.sampleMissingCount[i]; sampleGenotypedCount[i] += ls.sampleGenotypedCount[i]; sampleHetCount[i] += ls.sampleHetCount[i];
				sampleHomoRefCount[i] += ls.sampleHomoRefCount[i]; sampleHomoAltCount[i] += ls.sampleHomoAltCount[i]; sampleNonRefCount[i] += ls.sampleNonRefCount[i];
				sampleTsCount[i] += ls.sampleTsCount[i]; sampleTvCount[i] += ls.sampleTvCount[i]; sampleRareAlleleCount[i] += ls.sampleRareAlleleCount[i];
				sampleTotalDepth[i] += ls.sampleTotalDepth[i]; sampleDepthCount[i] += ls.sampleDepthCount[i];
			}

			for (int i = 0; i < 10; i++) siteMissingnessHistogram[i] += ls.siteMissingnessHistogram[i];
			for (int i = 0; i < 50; i++) mafHistogram[i] += ls.mafHistogram[i];
			for (int i = 0; i < 10; i++) ehHistogram[i] += ls.ehHistogram[i];
			for (int i = 0; i < 5; i++) hweChiSqHistogram[i] += ls.hweChiSqHistogram[i];

			ls.variantsPerChromosome.forEach((chrom, count) -> variantsPerChromosome.merge(chrom, count, Integer::sum));
			ls.binnedDensity.forEach((chrom, bins) -> {
				Map<Integer, Integer> globalBins = binnedDensity.computeIfAbsent(chrom, k -> new HashMap<>());
				bins.forEach((bin, count) -> globalBins.merge(bin, count, Integer::sum));
			});

			if (fstAccum != null && ls.fstAccum != null) {
				for (int pi = 0; pi < numPops; pi++) {
					for (int pj = 0; pj < numPops; pj++) {
						fstAccum[pi][pj][0] += ls.fstAccum[pi][pj][0];
						fstAccum[pi][pj][1] += ls.fstAccum[pi][pj][1];
						fstAccum[pi][pj][2] += ls.fstAccum[pi][pj][2];
					}
				}
			}
		}

		meanEH     = ehCount > 0 ? ehSum / ehCount : 0.0;
		meanChiSq  = numHweTested > 0 ? chiSqSum / numHweTested : 0.0;
		sampleObsHet = new double[numSamples];
		sampleFstat  = new double[numSamples];
		for (int i = 0; i < numSamples; i++) {
			double oh = sampleGenotypedCount[i] > 0 ? (double) sampleHetCount[i] / sampleGenotypedCount[i] : 0.0;
			sampleObsHet[i] = oh;
			sampleFstat[i]  = meanEH > 0 ? 1.0 - oh / meanEH : 0.0;
		}

		if (tajimaSiteCount > 0 && numSegSites > 1) {
			piHat = piPerSiteSum / tajimaSiteCount; 
			int nEff = (int) Math.round((double) totalAllelesForTajima / tajimaSiteCount);
			if (nEff >= 2) {
				int S = numSegSites;
				double a1 = 0, a2 = 0;
				for (int i = 1; i < nEff; i++) { a1 += 1.0 / i; a2 += 1.0 / ((double)i * i); }
				thetaW = S / a1;
				double b1 = (nEff + 1.0) / (3.0 * (nEff - 1));
				double b2 = 2.0 * (nEff * nEff + nEff + 3.0) / (9.0 * nEff * (nEff - 1));
				double c1 = b1 - 1.0 / a1; double c2 = b2 - (nEff + 2.0) / (a1 * nEff) + a2 / (a1 * a1);
				double e1 = c1 / a1; double e2 = c2 / (a1 * a1 + a2);
				double varD = e1 * S + e2 * S * (S - 1);
				if (varD > 0) tajimaD = (piPerSiteSum - thetaW) / Math.sqrt(varD);
			}
		}

		if (fstAccum != null) {
			pairwiseFst = new double[numPops][numPops];
			for (int pi = 0; pi < numPops; pi++) {
				for (int pj = pi + 1; pj < numPops; pj++) {
					double count = fstAccum[pi][pj][2];
					if (count > 0) {
						double meanHt = fstAccum[pi][pj][0] / count;
						double meanHs = fstAccum[pi][pj][1] / count;
						pairwiseFst[pi][pj] = meanHt > 0 ? (meanHt - meanHs) / meanHt : 0;
						pairwiseFst[pj][pi] = pairwiseFst[pi][pj];
					}
				}
			}
		}
	}

	private void processLine(String line, int numSamples, LocalStats ls, int[] samplePopIdx) {
		String[] cols = line.split("\t");
		if (cols.length < 9 + numSamples) return;
		String chrom = cols[0]; int pos = Integer.parseInt(cols[1]); String ref = cols[3]; String alt = cols[4];
		String[] altAllelesList = alt.equals(".") ? new String[0] : alt.split(",");
		int numAllelesAtSite = 1 + altAllelesList.length;

		boolean isSnp = ref.length() == 1 && (altAllelesList.length == 0 || altAllelesList[0].length() == 1);
		if (isSnp) ls.numSnps++; else ls.numIndels++;

		if (alt.equals(".")) ls.numMonomorphic++;
		else if (altAllelesList.length == 1) ls.numBiallelic++;
		else ls.numMultiallelic++;

		ls.variantsPerChromosome.merge(chrom, 1, Integer::sum);
		ls.binnedDensity.computeIfAbsent(chrom, k -> new HashMap<>()).merge(pos / DENSITY_BIN_SIZE, 1, Integer::sum);

		int siteMissing = 0, totalGenotypedAlleles = 0;
		int[] siteAlleleCounts = new int[numAllelesAtSite];
		int n00 = 0, n01 = 0, n11 = 0;

		for (int i = 0; i < numSamples; i++) {
			String[] gData = cols[9 + i].split(":");
			String gt = gData[0];
			if (gt.startsWith(".")) { ls.sampleMissingCount[i]++; siteMissing++; }
			else {
				ls.sampleGenotypedCount[i]++;
				String[] alleles = gt.split("[/|]");
				boolean isHet = false; String firstAllele = alleles[0];
				for (String a : alleles) {
					if (!a.equals(firstAllele)) isHet = true;
					int idx = Integer.parseInt(a);
					if (idx < numAllelesAtSite) { siteAlleleCounts[idx]++; totalGenotypedAlleles++; }
				}
				if (isHet) ls.sampleHetCount[i]++;
				else { if (firstAllele.equals("0")) ls.sampleHomoRefCount[i]++; else ls.sampleHomoAltCount[i]++; }
				if (!firstAllele.equals("0") || isHet) ls.sampleNonRefCount[i]++;
				if (altAllelesList.length == 1 && alleles.length == 2) {
					if (gt.equals("0/0") || gt.equals("0|0")) n00++; else if (gt.equals("1/1") || gt.equals("1|1")) n11++; else n01++;
				}
			}
		}

		ls.siteMissingnessHistogram[Math.min(9, (int) (10.0 * siteMissing / numSamples))]++;

		if (totalGenotypedAlleles > 0) {
			int minAc = Integer.MAX_VALUE;
			for (int c : siteAlleleCounts) if (c > 0 && c < minAc) minAc = c;
			double maf = (double) minAc / totalGenotypedAlleles;
			ls.mafHistogram[Math.min(49, (int) (100.0 * maf))]++;

			double eh = 1.0;
			for (int c : siteAlleleCounts) eh -= Math.pow((double) c / totalGenotypedAlleles, 2);
			ls.ehSum += eh; ls.ehCount++;
			ls.ehHistogram[Math.min(9, (int) (10.0 * eh))]++;

			if (altAllelesList.length == 1 && isSnp) {
				if (isTransition(ref, altAllelesList[0])) ls.numTransitions++; else ls.numTransversions++;
				if (numSamples >= 5) {
					ls.numHweTested++;
					double pVal = HweUtils.calculateHweFisherExactTest(n00, n01, n11);
					if (pVal < 0.05) ls.numHweViolations++;
					double p = (double) (2 * n00 + n01) / (2 * (n00 + n01 + n11));
					double exp01 = 2.0 * p * (1.0 - p) * (n00 + n01 + n11);
					double chi2 = (exp01 > 0) ? Math.pow(n01 - exp01, 2) / exp01 : 0;
					ls.chiSqSum += chi2;
					ls.hweChiSqHistogram[chi2 < 2 ? 0 : chi2 < 3.84 ? 1 : chi2 < 7 ? 2 : chi2 < 10 ? 3 : 4]++;
				}
			}

			if (maf > 0) {
				ls.numSegSites++; ls.totalAllelesForTajima += totalGenotypedAlleles; ls.tajimaSiteCount++;
				double p = (double) siteAlleleCounts[0] / totalGenotypedAlleles;
				ls.piPerSiteSum += 2.0 * p * (1.0 - p) * totalGenotypedAlleles / (totalGenotypedAlleles - 1.0);
			}

			if (samplePopIdx != null && altAllelesList.length == 1) {
				double pGlobal = (double) siteAlleleCounts[0] / totalGenotypedAlleles;
				double Ht = 2.0 * pGlobal * (1.0 - pGlobal);
				for (int pi = 0; pi < ls.numPops; pi++) {
					for (int pj = pi + 1; pj < ls.numPops; pj++) {
						double pI = getPopFreq(cols, numSamples, pi, samplePopIdx);
						double pJ = getPopFreq(cols, numSamples, pj, samplePopIdx);
						if (!Double.isNaN(pI) && !Double.isNaN(pJ)) {
							double Hs = (2.0 * pI * (1 - pI) + 2.0 * pJ * (1 - pJ)) / 2.0;
							ls.fstAccum[pi][pj][0] += Ht; ls.fstAccum[pi][pj][1] += Hs; ls.fstAccum[pi][pj][2]++;
						}
					}
				}
			}
		}
	}

	private double getPopFreq(String[] cols, int numSamples, int popIdx, int[] samplePopIdx) {
		int count0 = 0, total = 0;
		for (int i = 0; i < numSamples; i++) {
			if (samplePopIdx[i] == popIdx) {
				String gt = cols[9 + i].split(":")[0];
				if (!gt.startsWith(".")) {
					for (String a : gt.split("[/|]")) { if (a.equals("0")) count0++; total++; }
				}
			}
		}
		return total > 0 ? (double) count0 / total : Double.NaN;
	}

	private static class LocalStats {
		int numSnps=0, numIndels=0, numTransitions=0, numTransversions=0, numHweViolations=0, numHweTested=0, numMonomorphic=0, numBiallelic=0, numMultiallelic=0, numSegSites=0, tajimaSiteCount=0, ehCount=0;
		double ehSum=0, chiSqSum=0, piPerSiteSum=0;
		long totalAllelesForTajima=0;
		int[] sampleMissingCount, sampleGenotypedCount, sampleHetCount, sampleHomoRefCount, sampleHomoAltCount, sampleNonRefCount, sampleTsCount, sampleTvCount, sampleRareAlleleCount, sampleDepthCount;
		long[] sampleTotalDepth;
		int[] siteMissingnessHistogram = new int[10], mafHistogram = new int[50], ehHistogram = new int[10], hweChiSqHistogram = new int[5];
		Map<String, Integer> variantsPerChromosome = new HashMap<>();
		Map<String, Map<Integer, Integer>> binnedDensity = new HashMap<>();
		double[][][] fstAccum;
		int numPops;
		LocalStats(int numSamples, int numPops) {
			this.numPops = numPops;
			sampleMissingCount = new int[numSamples]; sampleGenotypedCount = new int[numSamples]; sampleHetCount = new int[numSamples];
			sampleHomoRefCount = new int[numSamples]; sampleHomoAltCount = new int[numSamples]; sampleNonRefCount = new int[numSamples];
			sampleTsCount = new int[numSamples]; sampleTvCount = new int[numSamples]; sampleRareAlleleCount = new int[numSamples];
			sampleDepthCount = new int[numSamples]; sampleTotalDepth = new long[numSamples];
			if (numPops > 0) fstAccum = new double[numPops][numPops][3];
		}
	}

	private boolean isTransition(String ref, String alt) {
		ref = ref.toUpperCase(); alt = alt.toUpperCase();
		return (ref.equals("A") && alt.equals("G")) || (ref.equals("G") && alt.equals("A")) ||
			   (ref.equals("C") && alt.equals("T")) || (ref.equals("T") && alt.equals("C"));
	}
}
