# BioCenicana: Sequential Analysis Pipeline for Polyploid Genomics

BioCenicana is a high-performance Java toolkit optimized for **Saccharum spp. (Sugarcane)** and other complex polyploids. It uses a line-by-line streaming engine to process massive VCF files with minimal memory footprint.

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
  -o integrated_report.tsv \
  --viz synteny_explorer.html
```
### Command Parameters Explained
*   `--gff1` / `--gff2`: **(Required)** The GFF files used as input for McScanX (usually a simplified 4-column format: Chr, GeneID, Start, End). These provide the core genomic coordinates for the structural visualization.
*   `--collinearity`: **(Required)** The standard `.collinearity` output file generated by McScanX or SynMap. This file contains the actual syntenic block groupings and gene-to-gene alignments.
*   `--annot1` / `--annot2`: **(Optional)** Functional annotation files. You can provide a standard **GFF3** (the tool will extract the `description` or `Note` from the 9th column) or a simple **TSV** file (`GeneID \t Description`). This powers the Gene Search and Functional Tooltips.
*   `--vcf`: **(Optional)** A standard Variant Call Format file containing population-level SNP data. The engine uses fuzzy matching to map these variants to the syntenic block boundaries, generating an "Evolutionary Heatmap" based on Marker Abundance (SNPs/Kbp).
*   `-o`: **(Optional)** Path to generate a consolidated Tab-Separated Values (TSV) report linking syntenic blocks, gene coordinates, and sequence data.
*   `--viz`: **(Optional)** Path to generate the standalone interactive HTML dashboard (Advanced Synteny Explorer).

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
    *   **Biological Legend**: The UI includes a mathematical description (0.00 to Max SNPs/Kbp) and biological context to facilitate interpretation and publication presentation.

#### 3. Intelligent Search and Statistics
*   **Functional Integration (`--annot1`, `--annot2`)**: Injects extracted biological functions (GO, PFAM) directly into interactive tooltips on mouse hover.
*   **Semantic Search**: Allows visual isolation of a block by searching for a specific Gene ID or function (e.g., "Kinase domain"). Matching blocks will turn yellow, dimming the rest of the genome.
*   **Focus Metrics**: Clicking any chromosome bar filters the visualization exclusively to that chromosome's blocks, displaying instant statistics: Total conserved genes, % of inversions, and the **Top 3 Biological Functions** in the region.
*   **Direct Export**: Instantly capture a high-resolution vector SVG ready for papers or Illustrator.

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
