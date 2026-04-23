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
Link individual SNP behavior to the population structure calculated in Step 4.
```bash
java -jar target/biocenicana-1.0.jar snp-explorer --matrix dosages_raw.tsv --pca my_population.pca.csv -p 10 -o audit.html
```
*   **Visual Check**: Open `audit.html` to see how dosages cluster for every SNP and how they map geographically on the PCA plot.

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
