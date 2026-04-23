# BioCenicana: Bioinformatics Toolkit for Polyploid Genomics

BioCenicana is a high-performance Java toolkit designed for the analysis of complex polyploid genomes, specifically optimized for **Saccharum spp. (Sugarcane)**. The toolkit focuses on memory efficiency through streaming data processing, allowing the analysis of massive VCF files that would otherwise exceed standard RAM limits.

The toolkit provides "Elite" grade analysis for population structure, kinship, linkage disequilibrium, and quality control.

---

## 1. Technical Features and Architecture

*   **Polyploid Awareness**: All algorithms are designed for polysomic inheritance and multiallelic systems.
*   **Streaming Engine**: Processes VCF files line-by-line to minimize memory footprint.
*   **Interactive Results**: Generates standalone HTML dashboards using Plotly.js for immediate visual auditing.
*   **Multi-Algorithm Clustering**: Integration of K-Means, DBSCAN, and GMM for robust population group detection.

---

## 2. Installation and Requirements

### System Requirements
*   **Java Runtime Environment (JRE)**: Version 11 or higher.
*   **Apache Maven**: For building from source.

### Building from Source
1. **Clone the repository**:
   ```bash
   git clone https://github.com/jhtrujillo/biocenicana.git
   cd biocenicana
   ```

2. **Build the project**:
   ```bash
   mvn clean package
   ```

3. **Executable JAR**:
   The resulting artifact will be located in:
   `target/biocenicana-1.0.jar`

---

## 3. Subcommand Reference

### A. vcf-stats
Generates a comprehensive statistical profile of a VCF file.

**Execution Command:**
```bash
java -jar biocenicana-1.0.jar vcf-stats -v input.vcf -o output_dir -p 10
```

**Key Parameters:**
*   `-v, --vcf`: Path to the input VCF file.
*   `-o, --output`: Base name for the output folder and results.
*   `-p, --ploidy`: Ploidy level of the samples (optional for stats).
*   `--popmap`: Optional file with sample population assignments for Fst.

**Detailed Results:**
*   **Console Summary**: Total counts of SNPs and individuals.
*   **dashboard.html**: Interactive charts showing MAF spectrum, depth distribution, and missingness.
*   **stats.tsv**: Comprehensive table with Transition/Transversion (Ts/Tv) ratios and HWE metrics for every site.

---

### B. vcf-filter
Filter variants based on quality metrics and select the most informative markers.

**Execution Command:**
```bash
java -jar biocenicana.jar vcf-filter -v input.vcf -o filtered.vcf --min-maf 0.05 --max-missing 0.2 --top-n 5000
```

**Filtering Logic:**
*   **MAF Filtering**: Removes SNPs with allele frequency below the threshold.
*   **Missingness Control**: Discards sites with high proportions of missing data.
*   **Top-N Selection**: Heuristically extracts the $N$ most informative SNPs based on higher Expected Heterozygosity (He).

---

### C. allele-dosage
Calculates allele dosages from raw read depths or genotype tags.

**Standard Dosage Matrix:**
```bash
java -jar biocenicana.jar allele-dosage -v input.vcf -p 10 > matrix.tsv
```

**Raw Frequency Export (For Clustering):**
```bash
java -jar biocenicana.jar allele-dosage -v input.vcf -p 10 --raw > raw_dosages.tsv
```

**Technical Note**: The `--raw` flag is essential for auditing SNP quality, as it exports unrounded continuous values (0.0 to 1.0) instead of discrete integers.

---

### D. pop-structure
Comprehensive population structure and kinship module.

**Execution Command:**
```bash
java -jar biocenicana.jar pop-structure --vcf input.vcf --ploidy 10 --pcs 10 --output my_pop
```

**Generated Assets:**
*   **my_pop.pca.csv**: Table containing PC coordinates, cluster assignments (K-Means, GMM, DBSCAN), and Ancestry Q-proportions.
*   **my_pop.kinship.csv**: Genomic Relationship Matrix (GRM) calculated via the VanRaden algorithm.
*   **my_pop.pca.html**: Standalone dashboard with 3D/2D PCA plots, Admixture barplots, and a UPGMA Dendrogram.

---

### E. snp-explorer
Integrated diagnostic tool for auditing SNP dosage quality vs. population structure.

**Execution Command:**
```bash
java -jar biocenicana.jar snp-explorer --matrix raw_dosages.tsv --pca my_pop.pca.csv --ploidy 10 -o audit.html
```

**Interactive Panels:**
1.  **Histogram View**: Displays empirical centroids and dosage distribution for the selected SNP.
2.  **PCA Map View**: Colors the whole population PCA map according to the dosage levels of the current SNP.

---

### F. ld (Linkage Disequilibrium)
Calculates pairwise $r^2$ between markers and estimates LD decay.

**Execution Command:**
```bash
java -jar biocenicana.jar ld -v input.vcf -o ld_report --max-dist 200000
```

---

### G. External Format Exporters
*   **gwaspoly-export**: Generates files compatible with the GWASpoly R package.
*   **joinmap**: Converts and fixes data for JoinMap linkage mapping.

---

## 4. Standard Integrated Workflow

To perform a complete "Elite" analysis of a polyploid population:

1. **Filtering**: Extract the 5000 most informative SNPs.
   ```bash
   java -jar biocenicana.jar vcf-filter -v raw.vcf -o filtered.vcf --top-n 5000
   ```

2. **Structure**: Map the population space.
   ```bash
   java -jar biocenicana.jar pop-structure -v filtered.vcf -p 10 -o pop_results
   ```

3. **Dosage**: Export raw allele frequencies.
   ```bash
   java -jar biocenicana.jar allele-dosage -v filtered.vcf -p 10 --raw > raw_mat.tsv
   ```

4. **Audit**: Verify SNP quality in the context of the PCA map.
   ```bash
   java -jar biocenicana.jar snp-explorer --matrix raw_mat.tsv --pca pop_results.pca.csv --ploidy 10 -o explorer.html
   ```

---

*This software is licensed under the MIT License. Developed for Advanced Genomic Breeding.*
