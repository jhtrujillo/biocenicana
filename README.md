# BioCenicana: Bioinformatics Toolkit for Polyploid Genomics

BioCenicana is a high-performance Java toolkit designed for the analysis of complex polyploid genomes, specifically optimized for **Saccharum spp. (Sugarcane)**. The toolkit focuses on memory efficiency through streaming data processing, allowing the analysis of massive VCF files that would otherwise exceed standard RAM limits.

---

## 🚀 Quick Start Tutorial: SNP Quality Audit

This tutorial walks you through the "Elite" workflow: filtering your data, calculating population structure, and auditing individual SNPs using the interactive dashboard.

### Step 1: Quality Filtering
Select the most informative markers to reduce noise.
```bash
java -jar target/biocenicana-1.0.jar vcf-filter -v raw_data.vcf -o filtered.vcf --min-maf 0.05 --top-n 5000
```

### Step 2: Population Structure (PCA)
Calculate the global structure backgrounds.
```bash
java -jar target/biocenicana-1.0.jar pop-structure -v filtered.vcf -p 10 -o my_pop
```

### Step 3: Extract Continuous Dosages
Export raw allele frequencies for cluster auditing.
```bash
java -jar target/biocenicana-1.0.jar allele-dosage -v filtered.vcf -p 10 --raw > dosages_raw.tsv
```

### Step 4: Integrated SNP Exploration
Launch the dashboard to audit specific SNPs.
```bash
java -jar target/biocenicana-1.0.jar snp-explorer --matrix dosages_raw.tsv --pca my_pop.pca.csv -p 10 -o audit.html
```

---

## 🛠 Subcommand Reference (Full List)

### 1. vcf-stats
Generates an interactive HTML dashboard and a TSV summary.
```bash
java -jar target/biocenicana-1.0.jar vcf-stats -v input.vcf -o stats_out -p 10
```
*   **Outputs**: `stats_out/stats_out.html` (visuals) and `stats_out/stats_out.tsv` (data).
*   **Metrics**: Ts/Tv ratio, MAF distribution, Depth per sample, Missingness.

### 2. vcf-filter
Advanced variant selection and cleaning.
```bash
java -jar target/biocenicana-1.0.jar vcf-filter -v in.vcf -o out.vcf --min-maf 0.05 --max-missing 0.2
```
*   `--top-n`: Heuristic selection of the $N$ SNPs with the highest heterozygosity.

### 3. allele-dosage
Computes polyploid dosage matrices (0 to Ploidy).
```bash
java -jar target/biocenicana-1.0.jar allele-dosage -v input.vcf -p 10 > matrix.tsv
```
*   `--raw`: Exports values as 0.0-1.0 frequencies (useful for clustering).

### 4. pop-structure
Population structure, kinship, and ancestry.
```bash
java -jar target/biocenicana-1.0.jar pop-structure -v input.vcf -p 10 -o output
```
*   **Results**: PCA coordinates (CSV), Kinship matrix (CSV), and Interactive Dashboard (HTML).

### 5. snp-explorer
Visual audit of SNP quality linked to PCA.
```bash
java -jar target/biocenicana-1.0.jar snp-explorer --matrix dosages.tsv --pca pca.csv -p 10 -o explorer.html
```

### 6. ld (Linkage Disequilibrium)
Pairwise $r^2$ calculation and decay visualization.
```bash
java -jar target/biocenicana-1.0.jar ld -v input.vcf -o ld_results --max-dist 200000
```

### 7. genetic-distance
Calculates pairwise genetic distance between all samples.
```bash
java -jar target/biocenicana-1.0.jar genetic-distance -v input.vcf -p 10 > dist_matrix.tsv
```

### 8. gwaspoly-export
Converts VCF data to the format required by the GWASpoly R package.
```bash
java -jar target/biocenicana-1.0.jar gwaspoly-export -v input.vcf -p 10 -o gwas_ready.csv
```

### 9. joinmap
Converts VCF data to JoinMap CP format for linkage mapping.
```bash
java -jar target/biocenicana-1.0.jar joinmap --input input.loc --output output.loc --fix
```

---

## 📦 Installation & Compilation

### Requirements
*   **Java JRE**: 11 or higher.
*   **Maven**: For building from source.

### Building
```bash
mvn clean package -DskipTests
# The executable JAR will be: target/biocenicana-1.0.jar
```

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
