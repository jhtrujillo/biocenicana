#!/usr/bin/env bash
# run_officinarum_pipeline.sh
# Integración automática de la Genómica Comparativa de Caña de Azúcar (CC-01-1940 vs S. officinarum LA-Purple)
set -e

echo "====================================================================="
echo "   INICIANDO PIPELINE DE INTEGRACIÓN CC-01-1940 vs S. OFFICINARUM     "
echo "====================================================================="

# Crear las carpetas destino
mkdir -p genomica_comparativa/1940_vs_officinarum/tables
mkdir -p genomica_comparativa/1940_vs_officinarum/plots

# 1. Compilación del proyecto BioJava
echo -e "\n[Paso 1/5] Compilando el proyecto BioJava con Maven..."
mvn clean package -DskipTests

# 2. Cálculo de Ka/Ks
if [ ! -f genomica_comparativa/1940_vs_officinarum/kaks_1940_vs_officinarum.tsv ]; then
  echo -e "\n[Paso 2/5] Calculando presiones evolutivas (Ka/Ks)..."
  java -jar target/biojava.jar kaks-calc \
    --collinearity benchmarks/genomica_comparativa/1940_vs_soff/1940_vs_soLA.collinearity \
    --cds1 benchmarks/genomas/1940/CC-01-1940.cds.fna \
    --cds2 benchmarks/genomas/s_officinarum_la/Saccharum_officinarum_LA-Purple.cds.fna \
    -o genomica_comparativa/1940_vs_officinarum/kaks_1940_vs_officinarum.tsv
else
  echo -e "\n[Paso 2/5] Archivo Ka/Ks ya existe. Saltando cálculo para ahorrar tiempo."
fi

# 3. Integración multi-ómica con comp-gen
echo -e "\n[Paso 3/5] Integrando GFFs, Colinealidad, VCF, Ka/Ks y longitudes..."
java -jar target/biojava.jar comp-gen \
  --gff1 benchmarks/genomas/s_officinarum_la/Saccharum_officinarum_LA-Purple.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_soff/1940_vs_soLA.collinearity \
  --cds1 benchmarks/genomas/s_officinarum_la/Saccharum_officinarum_LA-Purple.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --prot1 benchmarks/genomas/s_officinarum_la/Saccharum_officinarum_LA-Purple.protein.faa \
  --prot2 benchmarks/genomas/1940/CC-01-1940.protein.faa \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --kaks genomica_comparativa/1940_vs_officinarum/kaks_1940_vs_officinarum.tsv \
  --viz genomica_comparativa/1940_vs_officinarum/visor_sintenia.html \
  -o genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --name1 "S. officinarum" \
  --name2 "CC 1940" \
  --organism Saccharum

# 4. Análisis de genes relacionados con sacarosa
echo -e "\n[Paso 4/5] Analizando genes asociados al metabolismo de sacarosa..."
./.venv/bin/python scripts/check_sucrose_genes.py \
  --gff1 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --gff2 benchmarks/genomas/s_officinarum_la/Saccharum_officinarum_LA-Purple.gff3 \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --kaks genomica_comparativa/1940_vs_officinarum/kaks_1940_vs_officinarum.tsv \
  --name1 "CC 1940" \
  --name2 "S. officinarum"

# 5. Generación de Tablas y Gráficos Interactivos
echo -e "\n[Paso 5/5] Generando tablas y gráficos interactivos complementarios..."

# A. Orientación de bloques
./.venv/bin/python scripts/block_orientation.py \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --output genomica_comparativa/1940_vs_officinarum/tables/block_orientation.tsv

# B. Rangos de bloques
./.venv/bin/python scripts/compute_block_ranges.py \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --output genomica_comparativa/1940_vs_officinarum/tables/block_ranges.tsv

# C. Extraer SVs (vcf de referencia CC 1940)
./.venv/bin/python scripts/sv_parser.py \
  benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  genomica_comparativa/1940_vs_officinarum/tables/sv_regions.bed

# D. Intersección de SNPs de Azúcar
./.venv/bin/python scripts/sugar_snp_intersect.py \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --sugar-ids data/sugar_gene_ids.txt \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --output genomica_comparativa/1940_vs_officinarum/tables/sugar_snp_overlap.tsv

# E. Enriquecimiento GO/KEGG
./.venv/bin/python scripts/go_enrichment.py \
  --orient genomica_comparativa/1940_vs_officinarum/tables/block_orientation.tsv \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --output genomica_comparativa/1940_vs_officinarum/tables/go_enrichment.tsv \
  --organism sbicolor

# F. Intersección de bloques y SVs
./.venv/bin/python scripts/intersect_sv_blocks.py \
  --sv genomica_comparativa/1940_vs_officinarum/tables/sv_regions.bed \
  --orient genomica_comparativa/1940_vs_officinarum/tables/block_orientation.tsv \
  --output genomica_comparativa/1940_vs_officinarum/tables/sv_block_overlap.tsv

# G. Gráficos interactivos en HTML (Plotly)
./.venv/bin/python scripts/interactive_plots.py \
  --orient genomica_comparativa/1940_vs_officinarum/tables/block_orientation.tsv \
  --ranges genomica_comparativa/1940_vs_officinarum/tables/block_ranges.tsv \
  --report genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv \
  --sugar-ids data/sugar_gene_ids.txt \
  --snp-ovl genomica_comparativa/1940_vs_officinarum/tables/sugar_snp_overlap.tsv \
  --go-tsv genomica_comparativa/1940_vs_officinarum/tables/go_enrichment.tsv \
  --out-dir genomica_comparativa/1940_vs_officinarum/plots \
  --name1 "CC 1940" \
  --name2 "S. officinarum"

# H. Cuantificación de inversiones
./.venv/bin/python scripts/quantify_inversions.py \
  --ranges genomica_comparativa/1940_vs_officinarum/tables/block_ranges.tsv \
  --orient genomica_comparativa/1940_vs_officinarum/tables/block_orientation.tsv \
  --out-dir genomica_comparativa/1940_vs_officinarum/tables

echo -e "\n====================================================================="
echo "   ¡INTEGRACIÓN COMPLETADA CON ÉXITO!                               "
echo "   - Reporte tabular: genomica_comparativa/1940_vs_officinarum/reporte_comparativo.tsv"
echo "   - Visor HTML interactivo: genomica_comparativa/1940_vs_officinarum/visor_sintenia.html"
echo "====================================================================="
