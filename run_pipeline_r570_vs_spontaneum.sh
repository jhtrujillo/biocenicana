#!/usr/bin/env bash
# run_pipeline_r570_vs_spontaneum.sh
# Integración automática de la Genómica Comparativa de Caña de Azúcar (R570 vs Saccharum spontaneum AP85-441)
set -e

echo "====================================================================="
echo "   INICIANDO PIPELINE DE INTEGRACIÓN (R570 vs S. spontaneum AP)      "
echo "====================================================================="

# 1. Compilación del proyecto BioJava
echo -e "\n[Paso 1/5] Compilando el proyecto BioJava con Maven..."
mvn clean package -DskipTests

# Crear la carpeta destino
mkdir -p genomica_comparativa/r570_vs_spontaneum
mkdir -p genomica_comparativa/r570_vs_spontaneum/tables
mkdir -p genomica_comparativa/r570_vs_spontaneum/plots

# 2. Cálculo de Ka/Ks
if [ ! -f genomica_comparativa/r570_vs_spontaneum/kaks_r570_vs_spontaneum.tsv ]; then
  echo -e "\n[Paso 2/5] Calculando presiones evolutivas (Ka/Ks)..."
  java -jar target/biojava.jar kaks-calc \
    --collinearity dataset_genomico/comparativas/R570_vs_Spont_sim/R570_vs_Spont.collinearity \
    --cds1 dataset_genomico/genomas/R570_sim/R570.cds.fa \
    --cds2 dataset_genomico/genomas/Spont_sim/Spont.cds.fa \
    -o genomica_comparativa/r570_vs_spontaneum/kaks_r570_vs_spontaneum.tsv
else
  echo -e "\n[Paso 2/5] Archivo Ka/Ks ya existe. Saltando cálculo para ahorrar tiempo."
fi

# 3. Integración multi-ómica con comp-gen
echo -e "\n[Paso 3/5] Integrando GFFs, Colinealidad, VCF, Ka/Ks y longitudes..."
java -jar target/biojava.jar comp-gen \
  --gff1 dataset_genomico/genomas/R570_sim/R570.gff \
  --gff2 dataset_genomico/genomas/Spont_sim/Spont.gff \
  --collinearity dataset_genomico/comparativas/R570_vs_Spont_sim/R570_vs_Spont.collinearity \
  --cds1 dataset_genomico/genomas/R570_sim/R570.cds.fa \
  --cds2 dataset_genomico/genomas/Spont_sim/Spont.cds.fa \
  --prot1 dataset_genomico/genomas/R570_sim/R570.cds.fa \
  --prot2 dataset_genomico/genomas/Spont_sim/Spont.cds.fa \
  --vcf dataset_genomico/genomas/CC-01-1940/CC-01-1940_sim.vcf \
  --kaks genomica_comparativa/r570_vs_spontaneum/kaks_r570_vs_spontaneum.tsv \
  --viz genomica_comparativa/r570_vs_spontaneum/visor_sintenia.html \
  -o genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --name1 "R570" \
  --name2 "S. spontaneum" \
  --organism Saccharum

# 4. Análisis de genes relacionados con sacarosa
echo -e "\n[Paso 4/5] Buscando y cuantificando genes de metabolismo/transporte de azúcar..."
./.venv/bin/python scripts/check_sucrose_genes.py \
  --gff1 dataset_genomico/genomas/R570_sim/R570.gff \
  --gff2 dataset_genomico/genomas/Spont_sim/Spont.gff \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --kaks genomica_comparativa/r570_vs_spontaneum/kaks_r570_vs_spontaneum.tsv \
  --name1 "R570" \
  --name2 "S. spontaneum"

# 5. Generación de Tablas y Gráficos Interactivos
echo -e "\n[Paso 5/5] Generando tablas y gráficos interactivos complementarios..."

# A. Orientación de bloques
./.venv/bin/python scripts/block_orientation.py \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --output genomica_comparativa/r570_vs_spontaneum/tables/block_orientation.tsv

# B. Rangos de bloques
./.venv/bin/python scripts/compute_block_ranges.py \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --output genomica_comparativa/r570_vs_spontaneum/tables/block_ranges.tsv

# C. Extraer SVs
./.venv/bin/python scripts/sv_parser.py \
  dataset_genomico/genomas/CC-01-1940/CC-01-1940_sim.vcf \
  genomica_comparativa/r570_vs_spontaneum/tables/sv_regions.bed

# D. Intersección de SNPs de Azúcar
./.venv/bin/python scripts/sugar_snp_intersect.py \
  --vcf dataset_genomico/genomas/CC-01-1940/CC-01-1940_sim.vcf \
  --sugar-ids data/sugar_gene_ids.txt \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --output genomica_comparativa/r570_vs_spontaneum/tables/sugar_snp_overlap.tsv

# E. Enriquecimiento GO/KEGG
./.venv/bin/python scripts/go_enrichment.py \
  --orient genomica_comparativa/r570_vs_spontaneum/tables/block_orientation.tsv \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --output genomica_comparativa/r570_vs_spontaneum/tables/go_enrichment.tsv \
  --organism sbicolor

# F. Intersección de bloques y SVs
./.venv/bin/python scripts/intersect_sv_blocks.py \
  --sv genomica_comparativa/r570_vs_spontaneum/tables/sv_regions.bed \
  --orient genomica_comparativa/r570_vs_spontaneum/tables/block_orientation.tsv \
  --output genomica_comparativa/r570_vs_spontaneum/tables/sv_block_overlap.tsv

# G. Gráficos interactivos en HTML (Plotly)
./.venv/bin/python scripts/interactive_plots.py \
  --orient genomica_comparativa/r570_vs_spontaneum/tables/block_orientation.tsv \
  --ranges genomica_comparativa/r570_vs_spontaneum/tables/block_ranges.tsv \
  --report genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv \
  --sugar-ids data/sugar_gene_ids.txt \
  --snp-ovl genomica_comparativa/r570_vs_spontaneum/tables/sugar_snp_overlap.tsv \
  --go-tsv genomica_comparativa/r570_vs_spontaneum/tables/go_enrichment.tsv \
  --out-dir genomica_comparativa/r570_vs_spontaneum/plots \
  --name1 "R570" \
  --name2 "S. spontaneum"

# H. Cuantificación de inversiones
./.venv/bin/python scripts/quantify_inversions.py \
  --ranges genomica_comparativa/r570_vs_spontaneum/tables/block_ranges.tsv \
  --orient genomica_comparativa/r570_vs_spontaneum/tables/block_orientation.tsv \
  --out-dir genomica_comparativa/r570_vs_spontaneum/tables

echo -e "\n====================================================================="
echo "   ¡PIPELINE EJECUTADO CON ÉXITO!                                    "
echo "   - Reporte tabular: genomica_comparativa/r570_vs_spontaneum/reporte_comparativo.tsv"
echo "   - Visor HTML interactivo: genomica_comparativa/r570_vs_spontaneum/visor_sintenia.html"
echo "   - Gráficos interactivos en: genomica_comparativa/r570_vs_spontaneum/plots/"
echo "====================================================================="
