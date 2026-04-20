# 🧬 BioCenicana

[![Java Version](https://img.shields.io/badge/Java-11%2B-blue.svg)](https://www.oracle.com/java/)
[![Maven](https://img.shields.io/badge/Build-Maven-orange.svg)](https://maven.apache.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**BioCenicana** is a high-performance, streaming bioinformatics CLI toolkit meticulously designed for the complex genomics of polyploid crops, with a primary focus on Sugarcane (*Saccharum spp.*). It handles massive VCF files with memory efficiency, providing "Elite" grade population genetics and structure analysis.

---

## 🚀 Key Features

### 1. 🏆 Advanced Population Structure & Kinship Suite (`pop-structure`)
The most comprehensive population structure module for polyploids. It goes beyond simple PCA to provide a multi-layered diagnostic of your germplasm.

*   **Multi-Method Clustering:** Compare results from **K-Means** (centroid-based), **DBSCAN** (density-based for outlier detection), and **Gaussian Mixture Models (GMM)** (probabilistic/elliptical clustering).
*   **DAPC (Discriminant Analysis of Principal Components):** Maximizes genetic separation between groups, providing the clearest view of population differentiation.
*   **Ancestry & Admixture:** Interactive stacked barplots showing the genomic composition of every individual (Membership probabilities).
*   **Genomic Relationship Matrix (Kinship):** Implementation of the **VanRaden algorithm** for polyploids. Essential for controlling relatedness in Mixed-Model GWAS.
*   **Hierarchical Clustering (UPGMA):** Interactive dendrograms representing genetic kinship and lineage branching.
*   **Genetic Divergence:** Automatic calculation of **Global $F_{st}$** and pairwise distance matrices.

**Usage Example:**
```bash
java -jar biocenicana.jar pop-structure \
  --vcf input.vcf \
  --ploidy 10 \
  --pcs 10 \
  --output my_analysis
```

### 2. 🔍 Linkage Disequilibrium Analysis (`ld-analysis`)
Analyze how genomic regions are inherited together.
*   **LD Decay Plots:** Visualize how $r^2$ decreases over physical distance.
*   **Multi-threaded Engine:** High-speed LD calculation for large marker sets.
*   **Interactive Heatmaps:** Identify recombination blocks and haplotype structure.

**Usage Example:**
```bash
java -jar biocenicana.jar ld-analysis -v input.vcf -o ld_results --max-dist 100000
```

### 3. 🛡️ Polyploid-Aware VCF Filtering (`vcf-filter`)
Streaming filter that avoids memory bottlenecks while cleaning complex datasets.
*   **Dosage Estimation:** Automatically converts read depths (`AD`/`BSDP`) to accurate allele dosages for any ploidy level.
*   **Top-N Polymorphism Selection:** Uses a Min-Heap to extract only the $N$ most informative markers (based on Expected Heterozygosity).
*   **HWE & MAF Filtering:** Precise math for multiallelic and polysomic systems.

**Usage Example:**
```bash
java -jar biocenicana.jar vcf-filter -v raw.vcf -o clean.vcf -p 10 -m 0.05 -x 0.2 --top-n 5000
```

### 4. 📊 Population Genetics Statistics (`vcf-stats`)
Calculates highly precise metrics and generates a stunning interactive dashboard.
*   **MAF Spectrum:** High-resolution distribution of allele frequencies.
*   **Exact HWE:** Fisher's exact test for polyploids.
*   **Sample Quality Control:** TS/TV ratios, heterozygosity, and missingness per sample.

---

## 🎨 Interactive Dashboards
BioCenicana doesn't just output CSVs; it generates **Production-Ready HTML Dashboards** using `Plotly.js`. 

*   **3D PCA Exploration:** Rotate and zoom into your population space.
*   **Dynamic Coloring:** Toggle between different clustering methods in real-time.
*   **Hover Diagnostics:** Identify specific varieties directly on the charts.
*   **Zoomable Heatmaps:** Explore the Kinship matrix and Distance matrix with high resolution.

---

## 🛠️ Installation & Build

Requires **Java 11+** and **Maven**.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/biocenicana.git
    cd biocenicana
    ```

2.  **Compile with Maven:**
    ```bash
    mvn clean package
    ```

3.  **Run the tool:**
    ```bash
    java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar --help
    ```

---

## 📂 Output Formats
*   `.pca.csv`: Comprehensive table with PCs, Cluster IDs, LDs, and Ancestry Q-proportions.
*   `.kinship.csv`: Square Genomic Relationship Matrix (TASSEL/GAPIT/GWASpoly compatible).
*   `.pca.html`: Interactive diagnostic dashboard.
*   `.vcf`: Filtered and/or dosage-corrected VCF file.

---

## ⚖️ License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
*Developed for Advanced Genomic Breeding in Sugarcane.*
