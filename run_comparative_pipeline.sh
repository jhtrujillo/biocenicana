#!/usr/bin/env bash
# run_comparative_pipeline.sh
# Integración automática de la Genómica Comparativa de Caña de Azúcar (CC-01-1940 vs R570)
set -e

RUN_SV=false
for arg in "$@"; do
  if [ "$arg" == "--run-sv" ]; then
    RUN_SV=true
  fi
done

echo "====================================================================="
echo "   INICIANDO PIPELINE DE INTEGRACIÓN DE GENÓMICA COMPARATIVA         "
echo "====================================================================="

# 1. Compilación del proyecto BioJava
echo -e "\n[Paso 1/4] Compilando el proyecto BioJava con Maven..."
mvn clean package -DskipTests

# Crear la carpeta destino
mkdir -p genomica_comparativa/r570

# 2. Cálculo de Ka/Ks
if [ ! -f genomica_comparativa/r570/kaks_1940_vs_r570.tsv ]; then
  echo -e "\n[Paso 2/4] Calculando presiones evolutivas (Ka/Ks)..."
  java -jar target/biojava.jar kaks-calc \
    --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
    --cds1 benchmarks/genomas/1940/CC-01-1940.cds.fna \
    --cds2 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
    -o genomica_comparativa/r570/kaks_1940_vs_r570.tsv
else
  echo -e "\n[Paso 2/4] Archivo Ka/Ks ya existe en genomica_comparativa/r570/kaks_1940_vs_r570.tsv. Saltando cálculo para ahorrar tiempo."
fi

# 3. Integración multi-ómica con comp-gen
echo -e "\n[Paso 3/4] Integrando GFFs, Colinealidad, VCF, Ka/Ks y longitudes..."
java -jar target/biojava.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --cds1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --prot1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.protein.faa \
  --prot2 benchmarks/genomas/1940/CC-01-1940.protein.faa \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --kaks genomica_comparativa/r570/kaks_1940_vs_r570.tsv \
  --viz genomica_comparativa/r570/visor_sintenia.html \
  -o genomica_comparativa/r570/reporte_comparativo.tsv \
  --name1 "R570" \
  --name2 "CC 1940" \
  --organism Saccharum

# 4. Análisis de genes relacionados con sacarosa
./.venv/bin/python scripts/check_sucrose_genes.py \
  --gff1 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --gff2 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --kaks genomica_comparativa/r570/kaks_1940_vs_r570.tsv \
  --name1 "CC 1940" \
  --name2 "R570"

# 5. Generación de Tablas y Gráficos Interactivos
echo -e "\n[Paso 5/5] Generando tablas y gráficos interactivos complementarios..."

# A. Orientación de bloques
./.venv/bin/python scripts/block_orientation.py \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --output genomica_comparativa/r570/tables/block_orientation.tsv

# B. Rangos de bloques
./.venv/bin/python scripts/compute_block_ranges.py \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --output genomica_comparativa/r570/tables/block_ranges.tsv

# C. Extraer SVs (vcf de referencia CC 1940)
./.venv/bin/python scripts/sv_parser.py \
  benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  genomica_comparativa/r570/tables/sv_regions.bed

# D. Intersección de SNPs de Azúcar
./.venv/bin/python scripts/sugar_snp_intersect.py \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --sugar-ids data/sugar_gene_ids.txt \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --output genomica_comparativa/r570/tables/sugar_snp_overlap.tsv

# E. Enriquecimiento GO/KEGG
./.venv/bin/python scripts/go_enrichment.py \
  --orient genomica_comparativa/r570/tables/block_orientation.tsv \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --output genomica_comparativa/r570/tables/go_enrichment.tsv \
  --organism sbicolor

# F. Intersección de bloques y SVs
./.venv/bin/python scripts/intersect_sv_blocks.py \
  --sv genomica_comparativa/r570/tables/sv_regions.bed \
  --orient genomica_comparativa/r570/tables/block_orientation.tsv \
  --output genomica_comparativa/r570/tables/sv_block_overlap.tsv

# G. Gráficos interactivos en HTML (Plotly)
./.venv/bin/python scripts/interactive_plots.py \
  --orient genomica_comparativa/r570/tables/block_orientation.tsv \
  --ranges genomica_comparativa/r570/tables/block_ranges.tsv \
  --report genomica_comparativa/r570/reporte_comparativo.tsv \
  --sugar-ids data/sugar_gene_ids.txt \
  --snp-ovl genomica_comparativa/r570/tables/sugar_snp_overlap.tsv \
  --go-tsv genomica_comparativa/r570/tables/go_enrichment.tsv \
  --out-dir genomica_comparativa/r570/plots \
  --name1 "CC 1940" \
  --name2 "R570"

# H. Cuantificación de inversiones
./.venv/bin/python scripts/quantify_inversions.py \
  --ranges genomica_comparativa/r570/tables/block_ranges.tsv \
  --orient genomica_comparativa/r570/tables/block_orientation.tsv \
  --out-dir genomica_comparativa/r570/tables

echo -e "\n====================================================================="
echo "   ¡INTEGRACIÓN COMPLETADA CON ÉXITO!                               "
echo "   - Reporte tabular: genomica_comparativa/r570/reporte_comparativo.tsv"
echo "   - Visor HTML interactivo: genomica_comparativa/r570/visor_sintenia.html"
echo "====================================================================="

