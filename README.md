# BioCenicana

A high-performance, streaming bioinformatics CLI toolkit optimized for the complex genomics of Sugarcane (*Saccharum spp.*) and other polyploid crops.

## Features & Commands

BioCenicana provides several tools to analyze, filter, and process VCF files efficiently without memory bottlenecks.

### 1. `vcf-stats` (Advanced Statistics & HTML Dashboard)
Calculates highly precise population genetics statistics in a single pass. 
* **Multiallelic Expected Heterozygosity (EH):** Uses exact formulas (`1 - Σ(pi²)`) for accurate metrics in complex genomes.
* **Fisher Exact Test for HWE:** Mathematically exact p-values for Hardy-Weinberg Equilibrium, replacing standard Chi-square approximations.
* **Per-Sample Metrics:** Homozygous/Heterozygous counts, TS/TV ratios, rare allele counts.
* **Interactive HTML Dashboard:** Generates a stunning Plotly-based HTML report with high-resolution MAF spectra and HWE distribution charts.

**Usage:**
```bash
java -jar biocenicana.jar vcf-stats -v input.vcf -o output_dir
```

### 2. `vcf-filter` (Polyploid-Aware VCF Filtering)
A memory-efficient streaming filter that cleans your VCF files based on strict quality thresholds. Exceeds standard tools by supporting polysomic inheritance natively.

**Standard Filtering:**
```bash
java -jar biocenicana.jar vcf-filter -v input.vcf -o clean.vcf -m 0.05 -x 0.2 -e 0.05 -b
```
* `-m` / `--min-maf`: Minimum Minor Allele Frequency.
* `-x` / `--max-missingness`: Maximum allowed missing data (e.g. 0.2 = 20%).
* `-e` / `--min-hwe-pvalue`: Filter extreme HWE violations (e.g. 0.05).
* `-b` / `--biallelic-only`: Keeps only biallelic SNPs.

**Polyploid Mode (`-p` / `--ploidy`):**
Essential for sugarcane! If ploidy > 2, the tool ignores erroneous `0/1` diploid calls from GATK/NGSEP/FreeBayes and **automatically estimates allele dosages** using the `AD` or `BSDP` read depths. HWE and MAF are calculated accurately assuming polysomic inheritance.
```bash
java -jar biocenicana.jar vcf-filter -v cane.vcf -o clean.vcf -p 10 -e 0.05
```

**Top-N Polymorphic Selection (`-t` / `--top-n`):**
Perfect for designing targeted sequencing panels or SNP arrays. Uses an ultra-fast `PriorityQueue` (Min-Heap) to stream the VCF and extract only the top $N$ most polymorphic markers based on their Expected Heterozygosity (EH).
```bash
java -jar biocenicana.jar vcf-filter -v cane.vcf -o top_5000.vcf -p 10 -t 5000 --min-eh 0.4
```

### 3. Other Commands
* `allele-dosage`: Calculates allele dosages.
* `genetic-distance`: Computes genetic distance matrices.
* `joinmap`: Fixes/formats files for JoinMap.

## Compilation
Requires Java 11+ and Maven.
```bash
mvn clean package
```
*The executable fat JAR will be located at `target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar`.*
