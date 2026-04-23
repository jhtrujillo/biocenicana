# BioCenicana - Bioinformatics Toolkit for Polyploid Genomics

BioCenicana is a high-performance Java toolkit designed for the analysis of complex polyploid genomes, specifically optimized for Saccharum spp. (Sugarcane). The toolkit focuses on memory efficiency through streaming data processing, allowing the analysis of massive VCF files that would otherwise exceed standard RAM limits.

The toolkit provides "Elite" grade analysis for population structure, kinship, linkage disequilibrium, and quality control.

---

## Technical Features and Architecture

1. Polyploid Awareness: All algorithms are designed for polysomic inheritance and multiallelic systems.
2. Streaming Engine: Processes VCF files line-by-line to minimize memory footprint.
3. Interactive Results: Generates standalone HTML dashboards using Plotly.js for immediate visual auditing.
4. Multi-Algorithm Clustering: Integration of K-Means, DBSCAN, and GMM for robust population group detection.

---

## Installation and Requirements

Requirements:
- Java Runtime Environment (JRE) 11 or higher.
- Apache Maven (for building from source).

Building from source:
1. Clone the repository:
   git clone https://github.com/jhtrujillo/biocenicana.git
   cd biocenicana

2. Build the project:
   mvn clean package

3. The executable JAR will be located in:
   target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar

---

## Subcommand Reference

### 1. vcf-stats
Generates a comprehensive statistical profile of a VCF file.

Usage:
java -jar biocenicana.jar vcf-stats -v input.vcf -o output_directory --ploidy 10

Outputs:
- console summary: Basic counts of SNPs and samples.
- dashboard.html: Interactive charts showing MAF spectrum, depth distribution, and missingness.
- stats.tsv: Tab-separated table with metrics for every SNP.

Key Metrics:
- Major/Minor Allele Frequencies.
- Expected Heterozygosity (He).
- Missingness per site and per sample.
- Transition/Transversion (Ts/Tv) ratio.

---

### 2. vcf-filter
Filter variants based on quality metrics and select the most informative markers.

Usage:
java -jar biocenicana.jar vcf-filter -v input.vcf -o filtered.vcf --min-maf 0.05 --max-missing 0.2 --top-n 5000

Parameters:
- --min-maf: Minimum Minor Allele Frequency.
- --max-missing: Maximum proportion of missing data allowed per SNP.
- --top-n: Heuristically select the top N most informative SNPs (highest heterozygosity).

Results:
- A new VCF file containing only the variants that passed the specified filters.

---

### 3. allele-dosage
Calculates allele dosages (number of copies of the alternative allele) for each sample and SNP.

Usage (Standard):
java -jar biocenicana.jar allele-dosage -v input.vcf -p 10 > dosage_matrix.tsv

Usage (Raw Frequencies for Clustering):
java -jar biocenicana.jar allele-dosage -v input.vcf -p 10 --raw > raw_matrix.tsv

Description:
- Supports input from NGSEP, GATK, and FreeBayes.
- Uses read depth counts (AD/BSDP) to estimate continuous dosages.
- The --raw flag exports unrounded frequencies (0.0 to 1.0), which is essential for downstream empirical clustering.

---

### 4. pop-structure
Advanced population structure analysis using Principal Component Analysis and discriminatory methods.

Usage:
java -jar biocenicana.jar pop-structure --vcf input.vcf --ploidy 10 --pcs 10 --output my_pca

Outputs:
- my_pca.pca.csv: Table with sample names, PC coordinates, and cluster assignments.
- my_pca.kinship.csv: Square kinship matrix (VanRaden algorithm).
- my_pca.pca.html: Interactive 3D/2D dashboard for population exploration.

Analysis Components:
- PCA: Standard population projection.
- GMM and K-Means: Automatic detection of populations.
- DAPC: Maximization of between-group variance.
- Admixture Plot: Genomic ancestry proportions in a stacked bar chart.
- Dendrogram: UPGMA tree based on genetic distance.

---

### 5. snp-explorer
The integrated diagnostic tool for auditing SNP dosage quality and linkage to population structure.

Usage:
java -jar biocenicana.jar snp-explorer --matrix raw_matrix.tsv --pca my_pca.pca.csv --ploidy 10 --output snp_audit.html

Visualizations:
- Dosage Distribution Panel: Shows the 1D density and detected centroids for the selected SNP.
- Population PCA Panel: A full population PCA map where each sample is colored according to its dosage for the selected SNP.

Results:
- An explorer.html file that allows side-by-side comparison of local clustering vs global population structure.

---

### 6. ld (Linkage Disequilibrium)
Calculate pairwise r^2 between markers to understand recombination and decay.

Usage:
java -jar biocenicana.jar ld -v input.vcf -o ld_report --max-dist 200000

Outputs:
- ld_results.csv: Pairwise LD metrics.
- ld_decay.html: Visualization of LD reduction over physical distance.

---

### 7. genetic-distance
Calculates a pairwise genetic distance matrix (Manhattan/Euclidean) for large sets of individuals.

Usage:
java -jar biocenicana.jar genetic-distance -v input.vcf --ploidy 10 > distance_matrix.tsv

---

### 8. gwaspoly-export
Prepares the data for Genetic Association Studies using the GWASpoly package.

Usage:
java -jar biocenicana.jar gwaspoly-export -v input.vcf -p 10 -o gwas_input.csv

Description:
- Converts allele dosages into the ACGT/dosage format required for polyploid mixed models.

---

### 9. joinmap
Converts and fixes JoinMap CP format files for linkage mapping.

Usage:
java -jar biocenicana.jar joinmap --input data.loc --fix --output final.loc

---

## Standard Integrated Workflow Example

To perform a complete analysis of a new VCF dataset:

1. Clean the dataset:
   java -jar biocenicana.jar vcf-filter -v raw.vcf -o clean.vcf --min-maf 0.05 --top-n 10000

2. Calculate population structure:
   java -jar biocenicana.jar pop-structure -v clean.vcf -p 10 -o my_pop

3. Export raw frequencies for auditing:
   java -jar biocenicana.jar allele-dosage -v clean.vcf -p 10 --raw > dosages.tsv

4. Run the auditing explorer:
   java -jar biocenicana.jar snp-explorer --matrix dosages.tsv --pca my_pop.pca.csv --ploidy 10 -o audit.html

---

Developed for BioCenicana and the Sugarcane Genomics Community. This toolkit is licensed under the MIT License.
