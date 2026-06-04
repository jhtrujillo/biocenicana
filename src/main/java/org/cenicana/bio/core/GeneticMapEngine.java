package org.cenicana.bio.core;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * High-performance Genetic Mapping Engine.
 * Builds genetic linkage maps directly from population VCF files or segregation matrices.
 */
public class GeneticMapEngine {

    private final double minLod;
    private final double maxRecomb;
    private final String mappingFunction;
    private int ploidy = 10;
    private boolean singleDoseFilter = false;
    private double sdChi2PThreshold = 0.05;
    private int thinKb = 0;
    private boolean pseudoOnly = false;
    private int minLgMarkers = 1; // minimum markers per LG to include in output

    public void setThinKb(int thinKb) { this.thinKb = thinKb; }
    public void setPseudoOnly(boolean pseudoOnly) { this.pseudoOnly = pseudoOnly; }
    public void setMinLgMarkers(int minLgMarkers) { this.minLgMarkers = minLgMarkers; }

    private static final int MEMORY_WARN_THRESHOLD = 5000;
    private static final int MEMORY_HARD_LIMIT = 25000;

    public GeneticMapEngine(double minLod, double maxRecomb, String mappingFunction) {
        this.minLod = minLod;
        this.maxRecomb = maxRecomb;
        this.mappingFunction = mappingFunction != null ? mappingFunction.toLowerCase() : "kosambi";
    }

    public void setPloidy(int ploidy) {
        this.ploidy = ploidy;
    }

    /**
     * Activates the single-dose marker filter.
     * Only markers segregating 1:1 (simplex x nulliplex) are retained.
     * @param enabled  activate the filter
     * @param minChi2P minimum chi-square p-value for the 1:1 goodness-of-fit test
     */
    public void setSingleDoseFilter(boolean enabled, double minChi2P) {
        this.singleDoseFilter = enabled;
        this.sdChi2PThreshold = minChi2P;
    }

    public static class Marker {
        public String id;
        public String chr;
        public long pos;
        public double[] dosages; // NaN = missing

        public Marker(String id, String chr, long pos, int sampleCount) {
            this.id = id;
            this.chr = chr;
            this.pos = pos;
            this.dosages = new double[sampleCount];
        }
    }

    /**
     * Builds the genetic map from a VCF file and writes the results.
     */
    public void buildMap(String vcfPath, String mapOutputPath) throws IOException {
        // Validate output directory exists
        Path outPath = Paths.get(mapOutputPath);
        if (outPath.getParent() != null && !Files.exists(outPath.getParent())) {
            throw new IOException("Output directory does not exist: " + outPath.getParent());
        }

        System.out.println("[MapEngine] Reading variants from VCF: " + vcfPath);
        List<Marker> markers = parseVcf(vcfPath);
        if (markers.isEmpty()) {
            System.out.println("❌ Error: No markers loaded from VCF.");
            return;
        }
        int numMarkers = markers.size();
        System.out.println("[MapEngine] Loaded " + numMarkers + " markers across population.");

        if (numMarkers > MEMORY_HARD_LIMIT) {
            throw new IOException("Too many markers (" + numMarkers + "). Hard limit is " + MEMORY_HARD_LIMIT +
                    ". Apply MAF/missing filters to reduce the marker set before running genetic-map.");
        }
        if (numMarkers > MEMORY_WARN_THRESHOLD) {
            System.out.printf("[MapEngine] ⚠️  Large dataset (%d markers). Estimated memory: ~%.1f GB. This may be slow.%n",
                    numMarkers, 2.0 * numMarkers * numMarkers * 8.0 / 1e9);
        }

        // 1. Calculate Pairwise Recombination Frequencies and LOD scores (parallelized)
        double[][] rMatrix = new double[numMarkers][numMarkers];
        double[][] lodMatrix = new double[numMarkers][numMarkers];

        System.out.println("[MapEngine] Computing pairwise linkage (LOD and recombination frequencies)...");
        computePairwiseParallel(markers, rMatrix, lodMatrix);

        // 2. Group markers into Linkage Groups using Average-Linkage Clustering
        System.out.println("[MapEngine] Grouping markers into Linkage Groups (Chromosomes)...");
        int[] groups = partitionIntoLinkageGroups(numMarkers, lodMatrix, rMatrix);
        Map<Integer, List<Integer>> lgToMarkers = new HashMap<>();
        for (int i = 0; i < numMarkers; i++) {
            lgToMarkers.computeIfAbsent(groups[i], g -> new ArrayList<>()).add(i);
        }
        System.out.println("[MapEngine] Partition completed. Found " + lgToMarkers.size() + " Linkage Groups.");

        // 3. Order and position markers within each Linkage Group
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(mapOutputPath)))) {
            pw.println("Marker\tLinkageGroup\tPosition(cM)\tChr_Phys\tPos_Phys");

            // Sort LGs by size descending so larger, more informative groups get lower numbers
            List<Integer> sortedLgKeys = new ArrayList<>(lgToMarkers.keySet());
            sortedLgKeys.sort((a, b) -> Integer.compare(lgToMarkers.get(b).size(), lgToMarkers.get(a).size()));

            int lgIndex = 1;
            int discardedSmall = 0;
            for (int lgId : sortedLgKeys) {
                List<Integer> markerIndices = lgToMarkers.get(lgId);

                // Filter: skip LGs below minimum marker count
                if (markerIndices.size() < minLgMarkers) {
                    discardedSmall++;
                    continue;
                }

                if (markerIndices.size() < 2) {
                    Marker m = markers.get(markerIndices.get(0));
                    pw.printf(Locale.US, "%s\tLG%d\t0.00\t%s\t%d\n", m.id, lgIndex, m.chr, m.pos);
                    lgIndex++;
                    continue;
                }

                System.out.println("[MapEngine] Ordering LG" + lgIndex + " containing " + markerIndices.size() + " markers...");
                List<Integer> orderedPath = orderMarkersTSP(markerIndices, rMatrix);

                double currentPositionCm = 0.0;
                Marker firstMarker = markers.get(orderedPath.get(0));
                pw.printf(Locale.US, "%s\tLG%d\t0.00\t%s\t%d\n", firstMarker.id, lgIndex, firstMarker.chr, firstMarker.pos);

                for (int i = 1; i < orderedPath.size(); i++) {
                    int m1Idx = orderedPath.get(i - 1);
                    int m2Idx = orderedPath.get(i);
                    double r = rMatrix[m1Idx][m2Idx];
                    double distanceCm = convertRecombToCm(r);
                    currentPositionCm += distanceCm;

                    Marker m = markers.get(m2Idx);
                    pw.printf(Locale.US, "%s\tLG%d\t%.2f\t%s\t%d\n", m.id, lgIndex, currentPositionCm, m.chr, m.pos);
                }
                lgIndex++;
            }
            if (discardedSmall > 0)
                System.out.printf("[MapEngine] Filtered out %d LGs with fewer than %d markers.%n",
                        discardedSmall, minLgMarkers);
        }
        System.out.println("🎉 [MapEngine] Genetic linkage map successfully written to: " + mapOutputPath);
    }

    /**
     * Parses a VCF file to extract markers and dosage vectors.
     * Uses NaN for missing genotypes and resets dsIndex per variant line.
     * When singleDoseFilter is active, only simplex×nulliplex markers are retained.
     */
    private List<Marker> parseVcf(String vcfPath) throws IOException {
        List<Marker> markers = new ArrayList<>();
        String[] sampleNames = null;
        int discardedSingleDose = 0;
        int totalParsed = 0;

        try (BufferedReader br = Files.newBufferedReader(Paths.get(vcfPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("##")) continue;
                if (line.startsWith("#")) {
                    String[] cols = line.split("\t");
                    sampleNames = Arrays.copyOfRange(cols, 9, cols.length);
                    continue;
                }

                String[] cols = line.split("\t");
                if (cols.length < 10) continue;

                int numSamples = cols.length - 9;

                String chr = cols[0];
                long pos = Long.parseLong(cols[1]);
                String id = cols[2].equals(".") ? chr + "_" + pos : cols[2];

                // Resolve field indices fresh for every variant line.
                // Priority: ACN (Allele Copy Number, e.g. "9,1") > DS > GT
                int acnIndex = -1;
                int dsIndex  = -1;
                String[] formatFields = cols[8].split(":");
                for (int f = 0; f < formatFields.length; f++) {
                    if (formatFields[f].equalsIgnoreCase("ACN")) { acnIndex = f; break; }
                    if (formatFields[f].equalsIgnoreCase("DS"))  { dsIndex  = f; }
                }

                Marker marker = new Marker(id, chr, pos, numSamples);
                for (int i = 9; i < cols.length; i++) {
                    String sampleData = cols[i];
                    int sampleIdx = i - 9;

                    if (sampleData.startsWith(".")) {
                        marker.dosages[sampleIdx] = Double.NaN;
                        continue;
                    }

                    String[] values = sampleData.split(":");
                    double dosage;

                    if (acnIndex != -1 && values.length > acnIndex) {
                        // ACN format: "refCopies,altCopies" — use alt copy count as dosage
                        dosage = parseAcnDosage(values[acnIndex]);
                    } else if (dsIndex != -1 && values.length > dsIndex) {
                        try {
                            dosage = Double.parseDouble(values[dsIndex]);
                        } catch (NumberFormatException e) {
                            dosage = parseGenotypeDosage(values[0]);
                        }
                    } else {
                        dosage = parseGenotypeDosage(values[0]);
                    }
                    marker.dosages[sampleIdx] = dosage;
                }

                totalParsed++;

                // Single-dose filter: only keep simplex x nulliplex markers
                if (singleDoseFilter) {
                    SingleDoseResult sd = evaluateSingleDose(marker.dosages);
                    if (!sd.passes) {
                        discardedSingleDose++;
                        continue;
                    }
                }

                markers.add(marker);
            }
        }

        if (singleDoseFilter) {
            System.out.printf("[MapEngine] Single-dose filter: %d/%d markers retained (%.1f%%), %d discarded.%n",
                    markers.size(), totalParsed,
                    100.0 * markers.size() / Math.max(1, totalParsed),
                    discardedSingleDose);
        }

        // Filter: exclude contigs, keep only pseudochromosome markers
        if (pseudoOnly) {
            int before = markers.size();
            markers.removeIf(m -> m.chr.toLowerCase().contains("contig") ||
                                   m.chr.toLowerCase().contains("scaffold"));
            System.out.printf("[MapEngine] Pseudo-chromosome filter: %d/%d markers retained (contigs excluded).%n",
                    markers.size(), before);
        }

        // Physical thinning: keep 1 marker per thinKb kb window per chromosome
        if (thinKb > 0) {
            int before = markers.size();
            markers = thinByPosition(markers, thinKb * 1000L);
            System.out.printf("[MapEngine] Physical thinning (%d kb): %d/%d markers retained.%n",
                    thinKb, markers.size(), before);
        }

        return markers;
    }

    /**
     * Physical thinning: for each chromosome, sort markers by position and keep
     * only one per windowBp window (the first encountered after sorting).
     */
    private List<Marker> thinByPosition(List<Marker> markers, long windowBp) {
        // Group by chromosome
        Map<String, List<Marker>> byChr = new LinkedHashMap<>();
        for (Marker m : markers) byChr.computeIfAbsent(m.chr, k -> new ArrayList<>()).add(m);

        List<Marker> thinned = new ArrayList<>();
        for (List<Marker> chrMarkers : byChr.values()) {
            chrMarkers.sort(Comparator.comparingLong(m -> m.pos));
            long lastPos = -windowBp; // ensures first marker always passes
            for (Marker m : chrMarkers) {
                if (m.pos - lastPos >= windowBp) {
                    thinned.add(m);
                    lastPos = m.pos;
                }
            }
        }
        return thinned;
    }

    /**
     * Result of single-dose evaluation for one marker.
     */
    private static class SingleDoseResult {
        boolean passes;
        double chi2p;
        double freq;   // frequency of individuals with dosage >= 1
        double maxDose;
    }

    /**
     * Evaluates whether a marker qualifies as single-dose (simplex x nulliplex).
     *
     * Criteria:
     *  1. Maximum observed dosage == 1 (no individual carries more than one copy).
     *  2. Observed frequency of dosage-1 individuals is between 0.25 and 0.75
     *     (centered around the expected 0.5 for a 1:1 segregation).
     *  3. Chi-square goodness-of-fit for 1:1 (present vs absent) gives p >= sdChi2PThreshold,
     *     meaning the marker does NOT significantly deviate from 1:1.
     */
    private SingleDoseResult evaluateSingleDose(double[] dosages) {
        SingleDoseResult res = new SingleDoseResult();

        int nPresent = 0;   // dosage == 1
        int nAbsent  = 0;   // dosage == 0
        double maxObs = 0.0;
        int nValid = 0;

        for (double d : dosages) {
            if (Double.isNaN(d)) continue;
            nValid++;
            double rounded = Math.round(d);
            if (rounded > maxObs) maxObs = rounded;
            if (rounded == 1.0) nPresent++;
            else if (rounded == 0.0) nAbsent++;
            // dosage > 1 will make maxObs > 1 → criterion 1 fails
        }

        res.maxDose = maxObs;

        // Criterion 1: no multi-dose individuals
        if (maxObs > 1.0 || nValid < 5) {
            res.passes = false;
            return res;
        }

        // Criterion 2: frequency between 0.25 and 0.75
        res.freq = (double) nPresent / nValid;
        if (res.freq < 0.25 || res.freq > 0.75) {
            res.passes = false;
            return res;
        }

        // Criterion 3: chi-square test for 1:1 (only present vs absent considered)
        int nBinary = nPresent + nAbsent;
        if (nBinary < 5) { res.passes = false; return res; }

        double expected = nBinary / 2.0;
        double chi2 = Math.pow(nPresent - expected, 2) / expected
                    + Math.pow(nAbsent  - expected, 2) / expected;

        // chi2 CDF with df=1 using incomplete beta (reuse GwasMathUtils pattern inline)
        res.chi2p = 1.0 - chiSquareCDF(chi2, 1);
        res.passes = res.chi2p >= sdChi2PThreshold;
        return res;
    }

    /**
     * Chi-square CDF approximation using the incomplete gamma function.
     * For df=1: CDF(x) = erf(sqrt(x/2)).
     */
    private double chiSquareCDF(double x, int df) {
        if (x <= 0) return 0.0;
        if (df == 1) {
            return erf(Math.sqrt(x / 2.0));
        }
        // General case via regularized lower incomplete gamma (series expansion)
        double a = df / 2.0;
        double term = Math.exp(a * Math.log(x / 2.0) - (x / 2.0) - logGamma(a)) / a;
        double sum = term;
        for (int k = 1; k <= 200; k++) {
            term *= (x / 2.0) / (a + k);
            sum += term;
            if (term < sum * 1e-10) break;
        }
        return Math.min(1.0, sum);
    }

    private double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));
        double ans = 1.0 - t * Math.exp(-z * z - 1.26551223
                + t * (1.00002368 + t * (0.37409196 + t * (0.09678418
                + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398
                + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
        return z >= 0 ? ans : -ans;
    }

    private double logGamma(double x) {
        double[] c = {76.18009172947146, -86.50532032941677, 24.01409824083091,
                      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
        double y = x, tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        double ser = 1.000000000190015;
        for (double ci : c) ser += ci / ++y;
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    /**
     * Parses the ACN (Allele Copy Number) field to obtain alt allele dosage.
     * Format: "refCopies,altCopies" — e.g. "9,1" → dosage 1, "10,0" → dosage 0.
     * Returns NaN if the field is malformed or missing.
     */
    private double parseAcnDosage(String acn) {
        if (acn == null || acn.startsWith(".")) return Double.NaN;
        String[] parts = acn.split(",");
        if (parts.length < 2) return Double.NaN;
        try {
            return Double.parseDouble(parts[parts.length - 1].trim());
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    private double parseGenotypeDosage(String gt) {
        double d = 0.0;
        for (char c : gt.toCharArray()) {
            if (c == '1') d += 1.0;
        }
        return d;
    }

    /**
     * Parallel pairwise LOD and recombination frequency computation.
     * Spawns one task per row i, using all available CPU cores.
     */
    private void computePairwiseParallel(List<Marker> markers, double[][] rMatrix, double[][] lodMatrix)
            throws IOException {
        int numMarkers = markers.size();
        int threads = Runtime.getRuntime().availableProcessors();
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        AtomicInteger progress = new AtomicInteger(0);
        List<Future<?>> futures = new ArrayList<>();

        for (int i = 0; i < numMarkers; i++) {
            final int row = i;
            rMatrix[row][row] = 0.0;
            lodMatrix[row][row] = 100.0;
            futures.add(pool.submit(() -> {
                for (int j = row + 1; j < numMarkers; j++) {
                    double[] stats = calculatePairwiseLinkage(markers.get(row), markers.get(j));
                    rMatrix[row][j] = stats[0];
                    rMatrix[j][row] = stats[0];
                    lodMatrix[row][j] = stats[1];
                    lodMatrix[j][row] = stats[1];
                }
                int done = progress.incrementAndGet();
                if (done % Math.max(1, numMarkers / 20) == 0) {
                    System.out.printf("[MapEngine]   Pairwise progress: %d/%d rows (%.0f%%)%n",
                            done, numMarkers, 100.0 * done / numMarkers);
                }
            }));
        }

        pool.shutdown();
        try {
            for (Future<?> f : futures) f.get();
        } catch (InterruptedException | ExecutionException e) {
            pool.shutdownNow();
            throw new IOException("Parallel pairwise computation failed: " + e.getMessage(), e);
        }
    }

    /**
     * Computes recombination frequency (r) and LOD score between two markers.
     * LOD uses the standard genetics likelihood ratio formula, not a chi-square approximation.
     */
    private double[] calculatePairwiseLinkage(Marker m1, Marker m2) {
        int n = m1.dosages.length;
        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
        int count = 0;

        for (int k = 0; k < n; k++) {
            double v1 = m1.dosages[k];
            double v2 = m2.dosages[k];
            // FIX: use NaN check instead of sentinel -1.0
            if (Double.isNaN(v1) || Double.isNaN(v2)) continue;
            sumX += v1;
            sumY += v2;
            sumXY += v1 * v2;
            sumX2 += v1 * v1;
            sumY2 += v2 * v2;
            count++;
        }

        if (count < 5) {
            return new double[]{0.5, 0.0};
        }

        double mean1 = sumX / count;
        double mean2 = sumY / count;
        double num = sumXY - count * mean1 * mean2;
        double den1 = sumX2 - count * mean1 * mean1;
        double den2 = sumY2 - count * mean2 * mean2;

        if (den1 <= 0.0 || den2 <= 0.0) {
            return new double[]{0.5, 0.0};
        }

        double rCoef = num / Math.sqrt(den1 * den2);
        double rAbs = Math.abs(rCoef);
        double recomb = Math.max(0.0001, Math.min(0.4999, 0.5 * (1.0 - rAbs)));

        // FIX: standard genetics LOD = log10[ L(r) / L(0.5) ]
        // For biallelic marker pairs: LOD = n * [ r*log10(r) + (1-r)*log10(1-r) + log10(2) ]
        double lod = count * (recomb * Math.log10(recomb) + (1.0 - recomb) * Math.log10(1.0 - recomb) + Math.log10(2.0));
        lod = Math.max(0.0, lod);

        return new double[]{recomb, lod};
    }

    /**
     * Average-Linkage Hierarchical Clustering using an explicit member-list per group.
     * Each root tracks its members so average inter-group r/LOD is computed in O(|Gi|*|Gj|)
     * without scanning all n markers — avoids the O(n³) bottleneck on large datasets.
     */
    private int[] partitionIntoLinkageGroups(int numMarkers, double[][] lodMatrix, double[][] rMatrix) {
        int[] parent = new int[numMarkers];
        // Each root maps to its current member list
        Map<Integer, List<Integer>> groupMembers = new HashMap<>();
        for (int i = 0; i < numMarkers; i++) {
            parent[i] = i;
            List<Integer> members = new ArrayList<>();
            members.add(i);
            groupMembers.put(i, members);
        }

        // Collect only pairs that pass both thresholds, sorted by LOD descending
        System.out.println("[MapEngine] Building candidate link list...");
        List<int[]> links = new ArrayList<>();
        for (int i = 0; i < numMarkers; i++) {
            for (int j = i + 1; j < numMarkers; j++) {
                if (lodMatrix[i][j] >= minLod && rMatrix[i][j] <= maxRecomb) {
                    links.add(new int[]{i, j});
                }
            }
        }
        System.out.printf("[MapEngine] %d candidate links found. Merging groups...%n", links.size());
        links.sort((a, b) -> Double.compare(lodMatrix[b[0]][b[1]], lodMatrix[a[0]][a[1]]));

        int merged = 0;
        for (int[] link : links) {
            int ri = find(link[0], parent);
            int rj = find(link[1], parent);
            if (ri == rj) continue;

            List<Integer> groupI = groupMembers.get(ri);
            List<Integer> groupJ = groupMembers.get(rj);

            // Compute average r and LOD between the two groups
            double sumR = 0.0, sumLod = 0.0;
            int pairs = groupI.size() * groupJ.size();
            for (int a : groupI) {
                for (int b : groupJ) {
                    sumR   += rMatrix[a][b];
                    sumLod += lodMatrix[a][b];
                }
            }
            double avgR   = sumR   / pairs;
            double avgLod = sumLod / pairs;

            if (avgLod >= minLod && avgR <= maxRecomb) {
                // Merge smaller group into larger (union by size)
                int newRoot;
                if (groupI.size() >= groupJ.size()) {
                    parent[rj] = ri;
                    groupI.addAll(groupJ);
                    groupMembers.remove(rj);
                    newRoot = ri;
                } else {
                    parent[ri] = rj;
                    groupJ.addAll(groupI);
                    groupMembers.remove(ri);
                    newRoot = rj;
                }
                merged++;
                if (merged % 500 == 0) {
                    System.out.printf("[MapEngine]   Merged %d groups so far, %d groups remaining.%n",
                            merged, groupMembers.size());
                }
            }
        }

        System.out.printf("[MapEngine] Clustering complete: %d final groups.%n", groupMembers.size());
        int[] groups = new int[numMarkers];
        for (int i = 0; i < numMarkers; i++) {
            groups[i] = find(i, parent);
        }
        return groups;
    }

    private int find(int i, int[] parent) {
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i], parent);
    }

    private void union(int i, int j, int[] parent) {
        int rootI = find(i, parent);
        int rootJ = find(j, parent);
        if (rootI != rootJ) {
            parent[rootI] = rootJ;
        }
    }

    /**
     * TSP 2-Opt Heuristic to find optimal linear ordering of markers in a linkage group.
     * Uses delta-based swap evaluation (O(1) per swap) instead of recomputing the full path.
     */
    private List<Integer> orderMarkersTSP(List<Integer> markers, double[][] rMatrix) {
        // Nearest-neighbor initial path
        List<Integer> path = new ArrayList<>();
        Set<Integer> unvisited = new LinkedHashSet<>(markers);
        int current = markers.get(0);
        path.add(current);
        unvisited.remove(current);

        while (!unvisited.isEmpty()) {
            int next = -1;
            double minDist = Double.MAX_VALUE;
            for (int cand : unvisited) {
                if (rMatrix[current][cand] < minDist) {
                    minDist = rMatrix[current][cand];
                    next = cand;
                }
            }
            path.add(next);
            unvisited.remove(next);
            current = next;
        }

        // FIX: 2-opt with O(1) delta evaluation — only 2 edges change per swap
        boolean improved = true;
        int size = path.size();
        while (improved) {
            improved = false;
            for (int i = 0; i < size - 1; i++) {
                for (int j = i + 2; j < size; j++) {
                    // Current edges: (i, i+1) and (j, j+1 if exists)
                    int a = path.get(i);
                    int b = path.get(i + 1);
                    int c = path.get(j);
                    int d = (j + 1 < size) ? path.get(j + 1) : -1;

                    double currentCost = rMatrix[a][b] + (d != -1 ? rMatrix[c][d] : 0.0);
                    double swapCost    = rMatrix[a][c] + (d != -1 ? rMatrix[b][d] : 0.0);

                    if (swapCost < currentCost - 1e-10) {
                        // Reverse the segment [i+1 .. j] in place
                        int left = i + 1, right = j;
                        while (left < right) {
                            Collections.swap(path, left++, right--);
                        }
                        improved = true;
                    }
                }
            }
        }
        return path;
    }

    /**
     * Converts recombination frequency (r) into genetic distance (cM) using Kosambi or Haldane.
     */
    private double convertRecombToCm(double r) {
        r = Math.max(0.0001, Math.min(0.4999, r));
        if (mappingFunction.equalsIgnoreCase("haldane")) {
            return -50.0 * Math.log(1.0 - 2.0 * r);
        } else {
            // Kosambi
            return 25.0 * Math.log((1.0 + 2.0 * r) / (1.0 - 2.0 * r));
        }
    }
}
