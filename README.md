# BioJava: Sequential Analysis Pipeline for Polyploid Genomics

BioJava is a high-performance Java toolkit optimized for **Saccharum spp. (Sugarcane)** and other complex polyploids. It uses a line-by-line streaming engine to process massive VCF files with minimal memory footprint.

## Why BioJava?
*   **Streaming Engine**: Reads datasets line-by-line, decoupling RAM consumption from the number of variants. Run 60GB VCFs on a standard laptop.
*   **Polyploid-Aware Mathematics**: Discards discrete genotyping in favor of continuous allele dosages using Maximum Likelihood estimation.
*   **All-in-One Toolkit**: Go from raw VCF to PCA, Kinship, Phylogeny, LD decay, and Comparative Genomics in a single executable. No messy R scripts required.
*   **Interactive HTML Dashboards**: Every module exports a 100% offline HTML/D3.js interactive viewer. Democratize your results with non-bioinformatician breeders instantly.

## Quick Start (Try it in 2 minutes)
We provide a small simulated dataset so you can test BioJava immediately without needing your own data.

```bash
# 1. Clone the repository
git clone https://github.com/jhtrujillo/biojava.git
cd biojava

# 2. Compile the JAR file (Requires Java 11+ and Maven)
mvn clean package -DskipTests

# 3. Run the interactive PCA and Kinship module on the test VCF
java -jar target/biojava.jar pop-structure -v simulation_data/CC01_sim.vcf -p 10 -o test_output
```
*Open `test_output.pca.html` in your web browser to explore the 3D PCA plot!*

---

This manual follows the logical order of a standard bioinformatics pipeline.

---

## Step 1: Compilation and Setup
Before starting, ensure you have Java 11+ and Maven installed.
```bash
mvn clean package -DskipTests
# The executable JAR will be generated as: target/biojava.jar
```

---

## Step 2: Initial Dataset Diagnostics (`vcf-stats`)
Always start by understanding the raw state of your VCF.
```bash
java -jar target/biojava.jar vcf-stats -v raw_data.vcf -o initial_stats -p 10
```
*   **Result**: Creates an interactive dashboard (`initial_stats.html`) showing allele frequencies, depth distributions, and missingness.
*   **Use this to**: Decide your filtering thresholds (MAF and missingness).

---

## Step 3: Quality Control & Filtering (`vcf-filter`)
Clean your dataset to keep only high-quality, informative markers.
```bash
java -jar target/biojava.jar vcf-filter -v raw_data.vcf -o filtered.vcf --min-maf 0.05 --max-missing 0.2 --top-n 5000
```
*   **Result**: A new VCF (`filtered.vcf`) containing the top 5000 most heterozygous and complete SNPs.

---

## Step 4: Population Structure & Kinship (`pop-structure`)
Map the genetic space of your samples and calculate relationships.
```bash
java -jar target/biojava.jar pop-structure -v filtered.vcf -p 10 -o my_population
```
*   **Result**: Creates `my_population.pca.csv` (coordinates/clusters) and `my_population.kinship.csv` (VanRaden relationship matrix).
*   **Visualization**: Open `my_population.pca.html` to explore the population in 3D/2D and see ancestry barplots.

---

## Step 5: Dosage & Distance Matrix Extraction (`allele-dosage` & `genetic-distance`)
Export the finalized data for external statistical software.

**A. Allele Dosages (for GWAS/Mapping):**
```bash
java -jar target/biojava.jar allele-dosage -v filtered.vcf -p 10 --raw > dosages_raw.tsv
```

**B. Genetic Distance (for Phylogeny/Diversity):**
```bash
java -jar target/biojava.jar genetic-distance -v filtered.vcf -p 10 > matrix_distance.tsv
```

---

## Step 6: Interactive SNP Quality Audit (`snp-explorer`)
Audit individual SNP quality using **AD-Plots** (Reference vs Alternative depth scatter plots). By providing the VCF directly, the tool visualizes genotype clusters with high precision.
```bash
java -jar target/biojava.jar snp-explorer --vcf filtered.vcf --pca my_population.pca.csv --include list_of_snps.txt -p 10 -o audit.html
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
java -jar target/biojava.jar ld -v filtered.vcf -o ld_report --max-dist 200000
```
*   **Output**: An LD decay dashboard showing $r^2$ reduction over physical distance.

---

## Step 8: Exporting to Specialized Formats (`gwaspoly` & `joinmap`)
Finalize your analysis by connecting with other specialized tools.

**A. Export for R/GWASpoly:**
```bash
java -jar target/biojava.jar gwaspoly-export -v filtered.vcf -p 10 -o gwas_ready.csv
```

**B. Convert for JoinMap (Linkage Mapping):**
```bash
java -jar target/biojava.jar joinmap --input data.loc --output fixed.loc --fix
```

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*

---

## Step 9: Consolidating Batches (`vcf-merge`)
If you have data from different sequencing batches or chromosomes, use this to join them into a single master VCF. It automatically handles the union of samples and fills gaps with missing data.
```bash
java -jar target/biojava.jar vcf-merge -i batch1.vcf,batch2.vcf,batch3.vcf -o consolidated.vcf
```
*   **Intelligence**: Automatically detects samples in each file and creates a unified cross-table.

---

## Step 10: Comparative Genomics (`comp-gen`)
Integrate disparate genomic datasets to study synteny, structural variations, and functional conservation between genomes. This tool supports outputs from **McScanX** and **SynMap2 (CoGe)** and features an **Advanced Synteny Explorer** with integrated evolutionary analytics.

```bash
java -jar target/biojava.jar comp-gen \
  --gff1 genome1.gff \
  --gff2 genome2.gff \
  --collinearity results.collinearity \
  --annot1 annot1.tsv \
  --vcf population_variants.vcf \
  --organism Saccharum \
  --viz synteny_explorer.html
```

#### New Advanced Analytical Features
*   **Organism-Agnostic Rates**: Use `--organism` (Saccharum, Arabidopsis, Human, etc.) or `--subst-rate` to define custom substitution rates for divergence time estimation.
*   **Automated WGD Peak Detection**: Automatically identifies peaks in the $K_s$ distribution, estimates Millions of Years Ago (Mya) based on your organism's rate, and annotates them in the interactive histogram.
*   **Chromosome-Level Phylogeny**: Reconstructs a Neighbor-Joining (NJ) tree based on average $K_s$ distances between chromosomes, allowing you to visualize how subgenomes or individual chromosomes relate across species.
*   **BioJava AI Insights**: A built-in expert system that analyzes your current view (filters, genes, diversity) and generates a natural language report explaining the biological significance of the results.

#### Step 10.1: Native Ka/Ks Calculation (`kaks-calc`)
If you have your CDS sequences and a collinearity file, you can calculate the selection pressure (Ka/Ks ratios) natively using the Nei-Gojobori (1986) model:

```bash
java -jar target/biojava.jar kaks-calc \
  --collinearity simulation_data/R570_vs_CC01.collinearity \
  --cds1 simulation_data/R570.cds.fa \
  --cds2 simulation_data/CC01.cds.fa \
  -o real_kaks_results.tsv
```

Then, you can use the output file (`real_kaks_results.tsv`) as input for the `comp-gen` command with the `--kaks` flag to visualize it.

---

### Command Parameters Explained
*   `--gff1` / `--gff2`: **(Required)** The GFF files used as input for McScanX (usually a simplified 4-column format: Chr, GeneID, Start, End).
*   `--collinearity`: **(Required)** The standard `.collinearity` output file generated by McScanX or SynMap.
*   `--organism`: (Optional) Organism preset for substitution rate. Options: `Saccharum` (6.96e-9), `Arabidopsis` (7e-9), `Human` (1.2e-9), `Populus` (2.5e-9), `Drosophila` (1.5e-8). Default: `Saccharum`.
*   `--subst-rate`: (Optional) Custom substitution rate (subs/site/year). Overrides `--organism`.
*   `--annot1` / `--annot2`: (Optional) Functional annotation files (GFF3 or TSV).
*   `--vcf`: (Optional) Population variants for "Evolutionary Heatmap" (SNPs/Kbp).
*   `--kaks`: (Optional) A TSV file containing Ka/Ks substitution rates for gene pairs.
*   `-o`: (Optional) Path to generate a consolidated TSV report.
*   `--viz`: (Optional) Path to generate the standalone interactive HTML dashboard (Advanced Synteny Explorer).

### Advanced Synteny Explorer Features
By adding the `--viz` flag, the engine generates an interactive, high-performance HTML/D3.js dashboard with cutting-edge analytical tools:

#### 1. Multi-Modal Visualization (Chart Types)
*   **3D Ribbons**: Structurally connects Genome 1 and Genome 2 chromosomes using Bezier curves. Highlights macro-inversions in orange.
*   **Classic 2D Dotplot**: Projects sequences linearly (Genome 1 on the X-axis, Genome 2 on the Y-axis). Ideal for macroscopic overviews.
*   **Circular (Circos)**: Radial view to visualize inter-chromosomal relationships and translocations.

#### 2. Integrated Evolution (Color Modes)
*   **Evol. Heatmap**: Changes the entire chart palette to a thermal scale based on *marker abundance* (SNP Density/Kbp).
*   **Ka/Ks Selection**: Visualizes selection pressure on a **Green (Purifying) -> White (Neutral) -> Red (Positive)** scale.
*   **WGD Modal**: Interactive $K_s$ histogram with automated peak detection and real-time rate adjustment.

#### 3. Intelligent Analysis and AI
*   **AI Insights Assistant**: Click the "✨ AI Insights" button in the navbar to get a biological interpretation of your current selection, including functional enrichment and evolutionary stability analysis.
*   **Phylo Tree**: Explore the chromosome-level relationship tree with integrated biological summaries and distance metrics.
*   **Functional Integration**: Injects extracted biological functions (GO, PFAM) directly into interactive tooltips on mouse hover.
*   **Semantic Search**: Allows visual isolation of a block by searching for a specific Gene ID or function. Match results turn yellow.
*   **Direct Export**: Instantly capture a high-resolution vector SVG ready for papers.

---

## Step 11: SNP-Based Phylogenetic Tree (`snp-tree`)
Build an interactive **Neighbor-Joining phylogenetic tree** directly from a VCF file. The tool computes a pairwise genetic distance matrix and reconstructs the tree topology, then generates a standalone HTML viewer that works **100% offline** (no internet or external libraries required).

```bash
java -jar target/biojava.jar snp-tree \
  -v filtered.vcf \
  -p 10 \
  --output filogenia.nwk
```

*   **Output**: `filogenia.nwk` (Newick format) and `filogenia.html` (interactive viewer).

#### Optional: Color Nodes by PCA Clusters (`--pca`)
Combine the tree with population structure results from Step 4 (`pop-structure`) to color each leaf node by its cluster assignment:

```bash
java -jar target/biojava.jar snp-tree \
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

## Step 12: Advanced Genomic Annotation Dashboard (`annotate`)
The most advanced module in BioJava. It generates a professional-grade research dashboard that integrates Variant Effect Prediction (VEP), interactive synteny, population analysis, and laboratory-ready marker recommendations.

```bash
java -jar target/biojava.jar annotate \
  -v population.vcf \
  -g genome1.gff3 \
  -r genome1.fasta \
  --gff2 genome2.gff3 \
  -p proteins.faa \
  -c cds.fna \
  -w 5000 \
  -o annotation_suite.html
```

### Key Analytical Capabilities
*   **Hybrid Functional Search**: Combine a traditional functional catalog with a real-time autocomplete engine (GO, Pfam, InterPro) to isolate genes and markers instantly.
*   **Variant Effect Predictor (VEP)**: Automatically classifies mutations as **Missense**, **Synonymous**, **Stop-Gain**, or **Intronic** based on the primary GFF3 and CDS data.
*   **Interactive Synteny**: For every SNP, the tool visualizes the physical region in the primary genome and maps the orthologous block in a secondary genome (e.g., 1940 vs R570).
*   **KASP Marker Recommendation**: An expert system that scores SNPs (⭐ to ⭐⭐⭐) based on environment cleanliness (nearby SNPs) and GC content. Provides ready-to-order flanking sequences.
*   **Population Haplotype Matrix**: A global, color-coded grid visualizing genotypes across all samples and filtered markers to identify inheritance blocks.
*   **Sequence Explorer & BLAST**: Direct access to Protein and CDS sequences with integrated one-click BLAST functionality.

### Parameters
| Flag | Description |
|---|---|
| `-v` / `--vcf` | **(Required)** Input VCF file with samples |
| `-g` / `--gff` | **(Required)** Primary GFF3 file (Reference) |
| `-r` / `--ref-genome` | **(Required for KASP)** Reference Genome FASTA |
| `--gff2` | **(Optional)** Secondary GFF3 for Synteny |
| `-p` / `--protein` | **(Optional)** Protein FASTA file |
| `-c` / `--cds` | **(Optional)** CDS FASTA file |
| `-w` / `--window` | Search window in bp around markers (Default: 5000) |
| `-o` / `--output` | **(Required)** Output interactive HTML dashboard |

---

## Step 13: Polyploid GWAS Suite (`gwas`)
The most comprehensive module for genetic association in complex polyploids. It implements a high-performance **EMMAX (P3D)** linear mixed model engine capable of handling chromosome-specific kinship (LOCO), genomic windows (Haplotypes), and gene-gene interactions (Epistasis).

```bash
java -jar target/biojava.jar gwas \
  -v population.vcf \
  --pheno phenotypes.csv \
  --trait yield \
  --fixed location,year \
  --loco \
  --epistasis \
  -w 10 \
  -p 8 \
  -o gwas_dashboard.html
```

### Advanced Analytical Features
*   **EMMAX Engine**: Uses the P3D (Population Parameters Previously Determined) optimization to scan millions of markers in minutes.
*   **LOCO Method (`--loco`)**: *Leave-One-Chromosome-Out*. Calculates a unique kinship matrix for each chromosome, excluding the markers of the chromosome under analysis. This eliminates "proximal contamination" and significantly increases the power to detect true QTLs.
*   **Haplotype Blocks (`--window`)**: Instead of single SNPs, the engine groups $N$ markers into a sliding window and extracts the **First Principal Component (PC1)** via SVD to represent the local haplotype. Ideal for identifying regions under selection or blocks with high LD.
*   **Targeted Epistasis (`--epistasis`)**: Automatically identifies the top 5 "Lead" markers and scans the entire genome for $G \times G$ interactions (synergies). Identifies where the combined effect of two genes is greater than the sum of their parts.
*   **Categorical Fixed Effects (`--fixed`)**: Seamlessly integrates environmental variables (Locations, Replications, Years) into the model via automatic dummy encoding.
*   **Smart Contig Grouping**: In fragmented genomes (scaffolds/contigs), the engine automatically clusters small sequences to optimize LOCO performance without losing resolution.

### Interactive Discovery Dashboard
The output `.html` file provides a professional research console:
*   **Dual Manhattan View**: Toggle between a standard **Linear** view and a **Circular** (orbital) view for publication-quality genomic overviews.
*   **Interactive Selection Analysis**: Click on any significant marker to instantly see:
    *   **Boxplots**: Phenotype distribution by allele dosage.
    *   **Selection Guide**: Automatic recommendation (e.g., *"Select individuals with Dosage 4 to increase the trait"*).
    *   **Elite Candidates**: A dynamic table of individuals that carry the ideal genotype for that specific QTL.
*   **Genotype Reconstruction**: Displays the exact allelic configuration (e.g., `AAAT`, `GGCC`) based on the ploidy level.
*   **Epistasis Table**: A dedicated section for "Genetic Synergies" showing pairs of interacting genes.

### Parameters
| Flag | Description |
|---|---|
| `-v` / `--vcf` | **(Required)** Input VCF file |
| `--pheno` | **(Required)** Phenotype file (CSV/TSV). Must include a `Sample` column |
| `--trait` | **(Required)** Target trait name as found in the phenotype file |
| `-p` / `--ploidy` | **(Required)** Ploidy level (e.g., `4` for potato, `10` for sugarcane) |
| `--loco` | **(Optional)** Enable Leave-One-Chromosome-Out analysis (Recommended) |
| `-w` / `--window` | **(Optional)** Window size for Haplotype blocks (Default: 1, SNP-based) |
| `--epistasis` | **(Optional)** Run targeted GxG interaction scan for top hits |
| `--fixed` | **(Optional)** Comma-separated list of fixed effect columns (e.g., `env,block`) |
| `-k` / `--kinship` | **(Optional)** Path to pre-computed Kinship CSV |
| `--maf` | **(Optional)** Minimum Allele Frequency filter (Default: 0.05) |
| `-o` / `--output` | **(Optional)** Path for the interactive HTML dashboard |

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
