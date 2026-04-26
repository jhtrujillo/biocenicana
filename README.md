# BioCenicana: Sequential Analysis Pipeline for Polyploid Genomics

BioCenicana is a high-performance Java toolkit optimized for **Saccharum spp. (Sugarcane)** and other complex polyploids. It uses a line-by-line streaming engine to process massive VCF files with minimal memory footprint.

## Why BioCenicana?
*   **Streaming Engine**: Reads datasets line-by-line, decoupling RAM consumption from the number of variants. Run 60GB VCFs on a standard laptop.
*   **Polyploid-Aware Mathematics**: Discards discrete genotyping in favor of continuous allele dosages using Maximum Likelihood estimation.
*   **All-in-One Toolkit**: Go from raw VCF to PCA, Kinship, Phylogeny, LD decay, and Comparative Genomics in a single executable. No messy R scripts required.
*   **Interactive HTML Dashboards**: Every module exports a 100% offline HTML/D3.js interactive viewer. Democratize your results with non-bioinformatician breeders instantly.

## Quick Start (Try it in 2 minutes)
We provide a small simulated dataset so you can test BioCenicana immediately without needing your own data.

```bash
# 1. Clone the repository
git clone https://github.com/jhtrujillo/biocenicana.git
cd biocenicana

# 2. Compile the JAR file (Requires Java 11+ and Maven)
mvn clean package -DskipTests

# 3. Run the interactive PCA and Kinship module on the test VCF
java -jar target/biocenicana-1.0.jar pop-structure -v simulation_data/CC01_sim.vcf -p 10 -o test_output
```
*Open `test_output.pca.html` in your web browser to explore the 3D PCA plot!*

---

This manual follows the logical order of a standard bioinformatics pipeline.

---

## Step 1: Compilation and Setup
Before starting, ensure you have Java 11+ and Maven installed.
```bash
mvn clean package -DskipTests
# The executable JAR will be generated as: target/biocenicana-1.0.jar
```

---

## Step 2: Initial Dataset Diagnostics (`vcf-stats`)
Always start by understanding the raw state of your VCF.
```bash
java -jar target/biocenicana-1.0.jar vcf-stats -v raw_data.vcf -o initial_stats -p 10
```
*   **Result**: Creates an interactive dashboard (`initial_stats.html`) showing allele frequencies, depth distributions, and missingness.
*   **Use this to**: Decide your filtering thresholds (MAF and missingness).

---

## Step 3: Quality Control & Filtering (`vcf-filter`)
Clean your dataset to keep only high-quality, informative markers.
```bash
java -jar target/biocenicana-1.0.jar vcf-filter -v raw_data.vcf -o filtered.vcf --min-maf 0.05 --max-missing 0.2 --top-n 5000
```
*   **Result**: A new VCF (`filtered.vcf`) containing the top 5000 most heterozygous and complete SNPs.

---

## Step 4: Population Structure & Kinship (`pop-structure`)
Map the genetic space of your samples and calculate relationships.
```bash
java -jar target/biocenicana-1.0.jar pop-structure -v filtered.vcf -p 10 -o my_population
```
*   **Result**: Creates `my_population.pca.csv` (coordinates/clusters) and `my_population.kinship.csv` (VanRaden relationship matrix).
*   **Visualization**: Open `my_population.pca.html` to explore the population in 3D/2D and see ancestry barplots.

---

## Step 5: Dosage & Distance Matrix Extraction (`allele-dosage` & `genetic-distance`)
Export the finalized data for external statistical software.

**A. Allele Dosages (for GWAS/Mapping):**
```bash
java -jar target/biocenicana-1.0.jar allele-dosage -v filtered.vcf -p 10 --raw > dosages_raw.tsv
```

**B. Genetic Distance (for Phylogeny/Diversity):**
```bash
java -jar target/biocenicana-1.0.jar genetic-distance -v filtered.vcf -p 10 > matrix_distance.tsv
```

---

## Step 6: Interactive SNP Quality Audit (`snp-explorer`)
Audit individual SNP quality using **AD-Plots** (Reference vs Alternative depth scatter plots). By providing the VCF directly, the tool visualizes genotype clusters with high precision.
```bash
java -jar target/biocenicana-1.0.jar snp-explorer --vcf filtered.vcf --pca my_population.pca.csv --include list_of_snps.txt -p 10 -o audit.html
```
*   **Visual Check**: Open `audit.html`. This interactive dashboard allows you to:
    *   **Filter**: Use `--include` to focus only on specific SNPs of interest (one ID per line).
    *   **Histogram**: See the global distribution of allele frequencies.
    *   **AD-Plot**: See the Ref vs Alt depth scatter plot (Genotype clustering).
    *   **PCA View**: Map dosages onto the population structure from Step 4.

---

## Step 7: Linkage Disequilibrium Analysis (`ld`)
Study the recombination rates and LD decay.
```bash
java -jar target/biocenicana-1.0.jar ld -v filtered.vcf -o ld_report --max-dist 200000
```
*   **Output**: An LD decay dashboard showing $r^2$ reduction over physical distance.

---

## Step 8: Exporting to Specialized Formats (`gwaspoly` & `joinmap`)
Finalize your analysis by connecting with other specialized tools.

**A. Export for R/GWASpoly:**
```bash
java -jar target/biocenicana-1.0.jar gwaspoly-export -v filtered.vcf -p 10 -o gwas_ready.csv
```

**B. Convert for JoinMap (Linkage Mapping):**
```bash
java -jar target/biocenicana-1.0.jar joinmap --input data.loc --output fixed.loc --fix
```

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*

---

## Step 9: Consolidating Batches (`vcf-merge`)
If you have data from different sequencing batches or chromosomes, use this to join them into a single master VCF. It automatically handles the union of samples and fills gaps with missing data.
```bash
java -jar target/biocenicana-1.0.jar vcf-merge -i batch1.vcf,batch2.vcf,batch3.vcf -o consolidated.vcf
```
*   **Intelligence**: Automatically detects samples in each file and creates a unified cross-table.

---

## Step 10: Comparative Genomics (`comp-gen`)
Integrate disparate genomic datasets to study synteny, structural variations, and functional conservation between genomes. This tool supports outputs from **McScanX** and **SynMap2 (CoGe)** and features an Advanced Synteny Explorer.

```bash
java -jar target/biocenicana-1.0.jar comp-gen \
  --gff1 genome1.gff \
  --gff2 genome2.gff \
  --collinearity results.collinearity \
  --annot1 annot1.tsv \
  --annot2 annot2.gff3 \
  --vcf population_variants.vcf \
  --viz synteny_explorer.html
```

#### Step 10: Native Ka/Ks Calculation (`kaks-calc`)
If you have your CDS sequences and a collinearity file, you can calculate the selection pressure (Ka/Ks ratios) natively using the Nei-Gojobori (1986) model:

```bash
java -jar target/biocenicana-1.0.jar kaks-calc \
  --collinearity simulation_data/R570_vs_CC01.collinearity \
  --cds1 simulation_data/R570.cds.fa \
  --cds2 simulation_data/CC01.cds.fa \
  -o real_kaks_results.tsv
```

Then, you can use the output file (`real_kaks_results.tsv`) as input for the `comp-gen` command with the `--kaks` flag to visualize it.

---

### Command Parameters Explained
*   `--gff1` / `--gff2`: **(Required)** The GFF files used as input for McScanX (usually a simplified 4-column format: Chr, GeneID, Start, End). These provide the core genomic coordinates for the structural visualization.
*   `--collinearity`: **(Required)** The standard `.collinearity` output file generated by McScanX or SynMap. This file contains the actual syntenic block groupings and gene-to-gene alignments.
*   `--annot1` / `--annot2`: **(Optional)** Functional annotation files. You can provide a standard **GFF3** (the tool will extract the `description` or `Note` from the 9th column) or a simple **TSV** file (`GeneID \t Description`). This powers the Gene Search and Functional Tooltips.
*   `--vcf`: (Optional) A standard Variant Call Format file containing population-level SNP data. The engine uses fuzzy matching to map these variants to the syntenic block boundaries, generating an "Evolutionary Heatmap" based on Marker Abundance (SNPs/Kbp).
*   `--kaks`: (Optional) A TSV file containing Ka/Ks substitution rates for gene pairs (Format: Gene1 \t Gene2 \t Ka \t Ks \t Ka/Ks). Powers the "Selection Pressure" visualization mode.
*   `-o`: (Optional) Path to generate a consolidated Tab-Separated Values (TSV) report linking syntenic blocks, gene coordinates, and sequence data.
*   `--viz`: (Optional) Path to generate the standalone interactive HTML dashboard (Advanced Synteny Explorer).

### Advanced Synteny Explorer Features
By adding the `--viz` flag, the engine generates an interactive, high-performance HTML/D3.js dashboard with cutting-edge analytical tools:

#### 1. Multi-Modal Visualization (Chart Types)
*   **3D Ribbons**: Structurally connects Genome 1 and Genome 2 chromosomes using Bezier curves. Automatically detects block orientations and highlights chromosomal macro-inversions in orange.
*   **Classic 2D Dotplot**: Projects sequences linearly (Genome 1 on the X-axis, Genome 2 on the Y-axis). Ideal for macroscopic overviews, where descending diagonals mark collinear regions and ascending diagonals identify large-scale inversions. Switchable on the fly with a single click.

#### 2. Integrated Evolution (Color Modes)
*   **Structural Mode**: Colors blocks randomly (Ribbons) or highlights inverted blocks (Dotplot) to prioritize genome architecture analysis.
*   **Evol. Heatmap (Population Dynamics)**: Activated when a `--vcf` file is provided. Changes the entire chart palette to a thermal scale ("Turbo") based on *marker abundance* (SNP Density/Kbp):
    *   **Cold (Dark Blue)**: Ultra-conserved zones where no mutations were found in the population (Strong purifying selection).
    *   **Hot (Red/Orange)**: Polymorphism "Hotspots", revealing rapid evolution, recent mutations, or diversifying selection.
*   **Ka/Ks Selection (Selection Pressure)**: Activated when a `--kaks` file is provided. Visualizes the ratio of non-synonymous (Ka) to synonymous (Ks) substitutions using a diverging color scale:
    *   **Green (Ka/Ks < 1)**: Purifying selection. Indicates functional conservation and elimination of harmful mutations.
    *   **White (Ka/Ks ~ 1)**: Neutral evolution.
    *   **Red (Ka/Ks > 1)**: Positive selection. Indicates adaptive evolution or rapidly changing protein sequences.
    *   **Biological Legend**: The UI includes mathematical and biological context to facilitate interpretation and publication presentation.

#### 3. Intelligent Search and Statistics
*   **Functional Integration (`--annot1`, `--annot2`)**: Injects extracted biological functions (GO, PFAM) directly into interactive tooltips on mouse hover.
*   **Semantic Search**: Allows visual isolation of a block by searching for a specific Gene ID or function (e.g., "Kinase domain"). Matching blocks will turn yellow, dimming the rest of the genome.
*   **Focus Metrics**: Clicking any chromosome bar filters the visualization exclusively to that chromosome's blocks, displaying instant statistics: Total conserved genes, % of inversions, and the **Top 3 Biological Functions** in the region.
*   **Direct Export**: Instantly capture a high-resolution vector SVG ready for papers or Illustrator.

---

## Step 11: SNP-Based Phylogenetic Tree (`snp-tree`)
Build an interactive **Neighbor-Joining phylogenetic tree** directly from a VCF file. The tool computes a pairwise genetic distance matrix and reconstructs the tree topology, then generates a standalone HTML viewer that works **100% offline** (no internet or external libraries required).

```bash
java -jar target/biocenicana-1.0.jar snp-tree \
  -v filtered.vcf \
  -p 10 \
  --output filogenia.nwk
```

*   **Output**: `filogenia.nwk` (Newick format) and `filogenia.html` (interactive viewer).

#### Optional: Color Nodes by PCA Clusters (`--pca`)
Combine the tree with population structure results from Step 4 (`pop-structure`) to color each leaf node by its cluster assignment:

```bash
java -jar target/biocenicana-1.0.jar snp-tree \
  -v filtered.vcf \
  -p 10 \
  --output filogenia.nwk \
  --pca my_population.pca.csv
```

#### Interactive Viewer Features
Open the generated `.html` file in any browser to access:

*   **4 Layout Modes** (switchable on the fly):
    *   **Rectangular Cladogram**: Classical tree, all leaves at the same depth.
    *   **Rectangular Phylogram**: Branch lengths proportional to genetic distance.
    *   **Radial Cladogram**: Circular layout for topology overview.
    *   **Radial Phylogram**: Circular layout with real branch lengths.
*   **Genetic Distance Clustering**: A slider that dynamically cuts the tree at a distance threshold, automatically grouping and color-coding clades by genetic proximity.
*   **PCA Cluster Coloring** (requires `--pca`): Colors leaf nodes by their population structure cluster using any of the 3 methods available: **K-Means**, **DBSCAN**, or **GMM**.
*   **Sample Search**: Highlights matching nodes in red for quick identification.
*   **Interactive Tooltip**: Shows branch length and cumulative distance on hover.
*   **Pan / Zoom**: Mouse drag and scroll wheel for navigation.
*   **SVG Export**: Exports a high-resolution vector image for publications.

#### Parameters
| Flag | Description |
|---|---|
| `-v` / `--vcf` | **(Required)** Input VCF file |
| `-p` / `--ploidy` | **(Required)** Ploidy level (e.g., `10` for sugarcane) |
| `-o` / `--output` | **(Required)** Output Newick file path (`.nwk`) |
| `--pca` | **(Optional)** PCA CSV from `pop-structure` to color nodes by cluster |
| `-c` / `--caller` | Variant caller hint (`auto`, `gatk`, `freebayes`, `ngsep`). Default: `auto` |
| `-md` / `--min-depth` | Minimum genotype depth to trust a call. Default: `0` |
| `-t` / `--threads` | Number of CPU threads. Default: all available |

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
