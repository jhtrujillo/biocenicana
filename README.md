# BioCenicana: Bioinformatics Toolkit for Polyploid Genomics

BioCenicana is a high-performance Java toolkit designed for the analysis of complex polyploid genomes, specifically optimized for **Saccharum spp. (Sugarcane)**. The toolkit focuses on memory efficiency through streaming data processing, allowing the analysis of massive VCF files that would otherwise exceed standard RAM limits.

---

## 🚀 Quick Start Tutorial: SNP Quality Audit

This tutorial walks you through the "Elite" workflow: filtering your data, calculating population structure, and auditing individual SNPs using the interactive dashboard.

### Step 1: Quality Filtering
Start by selecting the most informative markers. This reduces noise and computational time.
```bash
java -jar target/biocenicana-1.0.jar vcf-filter -v raw_data.vcf -o filtered.vcf --min-maf 0.05 --top-n 5000
```
*   **Input**: Your raw VCF.
*   **Result**: `filtered.vcf` (The top 5000 most diverse and complete SNPs).

### Step 2: Population Structure (PCA)
Calculate the global structure of your population. This is the background for all downstream analyses.
```bash
java -jar target/biocenicana-1.0.jar pop-structure -v filtered.vcf -p 10 -o my_pop
```
*   **Result**: `my_pop.pca.csv` (Coordinates) and `my_pop.pca.html` (Interactive 3D plot).

### Step 3: Extract Continuous Dosages
Export the raw dosages (frequencies) for every sample. This is required for empirical cluster auditing.
```bash
java -jar target/biocenicana-1.0.jar allele-dosage -v filtered.vcf -p 10 --raw > dosages_raw.tsv
```
*   **Note**: Using `--raw` ensures we have the continuous distribution (0.0 to 1.0) instead of rounded integers.

### Step 4: Integrated SNP Exploration
Launch the dashboard to audit specific SNPs against the population structure.
```bash
java -jar target/biocenicana-1.0.jar snp-explorer --matrix dosages_raw.tsv --pca my_pop.pca.csv -p 10 -o audit.html
```
*   **Result**: `audit.html`. Open this file in your browser to see how each SNP's dosages group together and correlate with your PC1/PC2 map.

---

## Technical Features

*   **Polyploid Awareness**: Designed for polysomic inheritance and high ploidy levels (e.g., 8x, 10x, 12x).
*   **Streaming Engine**: Line-by-line processing for massive VCF files.
*   **Standalone Dashboards**: No server required; all visualizations are portable HTML files.

---

## Subcommand Reference

### vcf-stats
Generates a comprehensive statistical profile.
```bash
java -jar biocenicana-1.0.jar vcf-stats -v input.vcf -o output_dir -p 10
```
*   **Key Metrics**: MAF, Missingness per Sample/Site, Ts/Tv ratio, Heterozygosity.

### vcf-filter
Advanced variant selection.
*   `--min-maf`: Minimum Minor Allele Frequency.
*   `--max-missing`: Maximum missing data allowed per site.
*   `--top-n`: Selects the $N$ most informative markers.

### pop-structure
Population diagnostics and clustering.
*   **Clustering**: Automatic group detection using K-Means, GMM, and DBSCAN.
*   **Kinship**: VanRaden GRM matrix export.
*   **Ancestry**: Admixture barplots integrated in the dashboard.

### ld (Linkage Disequilibrium)
Pairwise $r^2$ and decay estimation.
```bash
java -jar biocenicana-1.0.jar ld -v input.vcf -o ld_report --max-dist 200000
```

---

## Installation

### Requirements
*   **Java JRE**: 11 or higher.
*   **Maven**: For building from source.

### Compilation
```bash
mvn clean package -DskipTests
# Result: target/biocenicana-1.0.jar
```

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
