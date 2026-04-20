# BioCenicana LD Decay Benchmarks

This folder contains Linkage Disequilibrium (LD) decay comparisons for different crops to demonstrate how genomic architecture affects the decay rate of $r^2$.

## Summary of Findings

| Crop | Breeding System | Decay Rate | Recommended Density |
| :--- | :--- | :--- | :--- |
| **Sugarcane** | Outcrossing / Polyploid | **Very Fast (~1 kb)** | High density (GBS/WGS) |
| **Maize** | Outcrossing / Diploid | **Moderate (~35 kb)** | Medium density (SNP Chips) |
| **Rice** | Selfing / Diploid | **Slow (~260 kb)** | Low density (Few thousand SNPs) |

## Folders

- `maize/`: Synthetic maize VCF and interactive dashboard.
- `rice/`: Synthetic rice VCF and interactive dashboard.
- `sugarcane/`: Placeholder for sugarcane results (Real data).

## How to reproduce
Use the BioCenicana `ld` command with appropriate window sizes:
- Sugarcane: `-w 50000 --bin-size 1000`
- Maize: `-w 200000 --bin-size 5000`
- Rice: `-w 1000000 --bin-size 20000`
