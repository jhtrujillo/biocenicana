# BioCenicana

BioCenicana is a high-performance bioinformatics CLI toolkit designed for the genomics of polyploid crops, with a primary focus on Sugarcane (Saccharum spp.). It handles massive VCF files with memory efficiency, providing advanced population genetics and structure analysis.

---

## Tool Overview and Functions

BioCenicana is organized into several specialized subcommands. Below is the complete list of available functions:

### 1. Population Structure Analysis (pop-structure)
Performs Principal Component Analysis (PCA) and detailed clustering.
- Multiple Clustering Methods: K-Means, DBSCAN, and Gaussian Mixture Models (GMM).
- DAPC (Discriminant Analysis of Principal Components): Enhances genetic separation between groups.
- Ancestry Proportions: Calculates and visualizes genomic composition (Admixture).
- Genomic Relationship Matrix (Kinship): Uses the VanRaden algorithm optimized for polyploids.
- Hierarchical Clustering: Generates UPGMA dendrograms.
- Global Fst: Calculates overall population differentiation.

### 2. SNP Quality and Group Explorer (snp-explorer)
An interactive visualization tool to audit the quality of allele dosages.
- Dosage Histograms: Visualize how dosages group across the population for each SNP.
- PCA Mapping: Dynamically color a population-wide PCA map using the dosage of a specific SNP.
- Performance: Generates a lightweight, interactive HTML report with Plotly.js.

### 3. Allele Dosage Calculation (allele-dosage)
Converts raw VCF data (AD, GT, BSDP) into discrete or continuous dosages.
- Raw Frequency Mode: Option to export unrounded frequencies (0.0 to 1.0) for clustering analysis.
- Flexible Parsing: Supports GATK, FreeBayes, NGSEP, and other standard VCF callers.

### 4. VCF Filtering (vcf-filter)
Advanced filtering logic designed for complex polyploid genomes.
- Informative SNPs: Selects the Top-N SNPs based on Expected Heterozygosity.
- Quality Filters: MAF (Minor Allele Frequency) and maximum missingness.
- Streaming Engine: Processes files of any size without high memory consumption.

### 5. Statistics and Quality Control (vcf-stats)
Generates summary statistics to evaluate the overall quality of the dataset.
- MAF Spectrum and HWE (Hardy-Weinberg Equilibrium).
- Sample Metrics: Heterozygosity, TS/TV ratio, and individual missingness.

### 6. Linkage Disequilibrium Analysis (ld-analysis)
Calculates LD ($r^2$) between markers and estimates LD decay over physical distance.

### 7. Genetic Distance (genetic-distance)
Computes pairwise genetic distance matrices between all individuals in the VCF.

### 8. External Format Export
- gwas-poly-export: Formats data for use in the GWASpoly R package.
- joinmap-export: Generates files compatible with JoinMap for genetic linkage mapping.

---

## Detailed Workflow: SNP + PCA Explorer

To use the new integrated visualization that shows dosages on top of the population structure, follow these steps:

### Step 1: Export Raw Dosages
Generate a tab-separated matrix of unrounded allele frequencies (0.0 - 1.0).
```bash
java -jar biocenicana.jar allele-dosage -v input.vcf -p 10 --raw > raw_matrix.tsv
```

### Step 2: Calculate Population Structure
Run the PCA analysis to generate the coordinate system for the population.
```bash
java -jar biocenicana.jar pop-structure -v input.vcf -p 10 -o results
```
This generates `results.pca.csv`, which contains the PC1 and PC2 coordinates.

### Step 3: Generate the Integrated Dashboard
Link the raw frequencies with the PCA coordinates to create the interactive explorer.
```bash
java -jar biocenicana.jar snp-explorer \
  --matrix raw_matrix.tsv \
  --pca results.pca.csv \
  --ploidy 10 \
  --output explorer.html
```

---

## Installation and Build

Requirements: Java 11 or higher and Maven.

1. Clone the repository:
```bash
git clone https://github.com/jhtrujillo/biocenicana.git
cd biocenicana
```

2. Compile:
```bash
mvn clean package
```

3. Run:
```bash
java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar --help
```

---

## License

This project is licensed under the MIT License. Developed for Advanced Genomic Breeding in Sugarcane.
