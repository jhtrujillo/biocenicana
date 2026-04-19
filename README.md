# Biocenicana

Bioinformatics tools for population genomics and molecular fingerprinting studies in polyploid organisms (such as sugarcane). Built with Java, Maven, and [PicoCLI](https://picocli.info/).

## Key Features

*   **`allele-dosage`**: Rapidly generates allele dosage matrices from massive VCF files, bypassing memory overhead and avoiding `OutOfMemoryError` exceptions. It features built-in support to auto-detect VCFs generated from NGSEP, GATK, and FreeBayes, allowing dynamic, on-the-fly imputation models.
*   **Optimized VCF Parsing**: Implements a highly efficient stream-based reader (`VcfFastReader.java`) to guarantee extremely low RAM consumption, which is ideal for processing highly repetitive and polyploid genomes.

## Requirements

*   **Java 11** or higher.
*   **Apache Maven 3.6+** (for building the project).

*(Note: The core `NGSEPcore` library dependency is included locally in the `/lib` folder of this repository).*

## Build Instructions

To build the project and generate the executable JAR file containing all dependencies, run the following command in the root directory:

```bash
mvn clean package -DskipTests
```

This will automatically create a standalone packaged jar in the `target/` folder named `biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar`.

## Usage

The packaged jar acts as a central executable toolkit containing multiple CLI commands. We recommend running it directly via `java -jar`.

### View general help menu
```bash
java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar --help
```

### Allele Dosage Matrix Calculation (`allele-dosage`)

Computes a matrix of genotypes converted to allele dosages for a given ploidy, extracting depths directly from the VCF FORMAT fields.

```bash
java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar allele-dosage \
  --vcf path/to/your_file.vcf \
  --ploidy 10 \
  --impute bsdp-mode
```

#### Available Options:
*   `-p, --ploidy`: Sets the ploidy level of the targeted organism (e.g., `10` for sugarcane). Defaults to `2`.
*   `-i, --impute`: Missing value imputation method.
    *   `bsdp`: Leaves missing/uncalled genotypes as `-1`.
    *   `mode`: Imputes missing genotypes with the mode dosage of the SNP across all samples.
    *   `mean`: Imputes missing genotypes with the mean dosage of the SNP across all samples.
    *   `bsdp-mode` / `bsdp-mean`: Dynamically extracts available counts first, and then imputes only the remaining `-1` values with the mode or mean.
*   `-c, --caller`: Define which Variant Calling software generated the VCF. Options: `auto`, `ngsep`, `gatk`, `freebayes`. If left empty, it defaults to `auto`.

## Code Structure

This repository was recently modernized to abandon monolithic scripts in favor of a clean, class-based architecture:
*   `src/main/java/org/cenicana/bio/cli/`: Contains the individual CLI endpoint commands managed by PicoCLI.
*   `src/main/java/org/cenicana/bio/io/`: Contains efficient custom parsers (like `VcfFastReader`) focused on I/O optimization.
*   `src/main/java/org/cenicana/bio/`: Legacy analytical classes (e.g., `GeneDosis.java`) refactored into modern data streams.
*   `lib/`: Contains critical external dependencies not found in Maven Central (e.g., `NGSEPcore`).

---
Project managed and maintained by [jhtrujillo](https://github.com/jhtrujillo) and the Cenicaña Bioinformatics Research Team.
