#!/usr/bin/env bash
# =============================================================================
# Pipeline: Mapa genético + ubicación de genes de sacarosa
# Actividad 4.1.5 — Asociación de genes a marcadores/regiones de interés
#
# Uso:
#   bash run_genetic_map_pipeline.sh
# =============================================================================

set -e

echo "============================================================"
echo "  PIPELINE: MAPA GENÉTICO + GENES DE SACAROSA              "
echo "============================================================"

# ── Rutas ────────────────────────────────────────────────────────────────────
VCF="benchmarks/vcfs/mapa_genetico/AllSamples_variants_geneticmap_pseudochromosomes_standarfilters_minInd_85_single_dosage.vcf"
GFF="benchmarks/genomas/1940/CC-01-1940.gff3"
GENES="data/sugar_gene_ids.txt"   # lista de genes de sacarosa (uno por línea)
OUT="genomica_comparativa/mapa_genetico"
PYTHON=".venv/bin/python"

# ── Parámetros del mapa ───────────────────────────────────────────────────────
PLOIDY=10
LOD=3.0
MAX_R=0.35
THIN_KB=500
MIN_LG=3
MAPPING_FUNC="kosambi"

mkdir -p "$OUT"

# ── Paso 1: Compilar BioJava ──────────────────────────────────────────────────
echo ""
echo "[Paso 1/3] Compilando BioJava..."
mvn clean package -DskipTests -q
echo "  ✓ target/biojava.jar listo"

# ── Paso 2: Construir el mapa genético ───────────────────────────────────────
echo ""
echo "[Paso 2/3] Construyendo mapa genético..."

java -Xmx8g -jar target/biojava.jar genetic-map \
  -i "$VCF" \
  -p $PLOIDY \
  -o "$OUT/mapa_biparental.map" \
  --pseudo-only \
  --thin-kb $THIN_KB \
  --min-lg-markers $MIN_LG \
  --lod $LOD \
  --max-r $MAX_R \
  --mapping-function $MAPPING_FUNC \
  --viz "$OUT/visor_mapa.html"

echo "  ✓ Mapa: $OUT/mapa_biparental.map"
echo "  ✓ Visor: $OUT/visor_mapa.html"

# ── Paso 3: Ubicar genes de sacarosa en el mapa ───────────────────────────────
echo ""
echo "[Paso 3/3] Ubicando genes de sacarosa en el mapa..."

$PYTHON scripts/map_genes_to_map.py \
  --map   "$OUT/mapa_biparental.map" \
  --gff   "$GFF" \
  --genes "$GENES" \
  --output "$OUT/genes_en_mapa.tsv"

echo ""
echo "============================================================"
echo "  RESULTADOS                                                "
echo "============================================================"
echo "  Mapa genético:      $OUT/mapa_biparental.map"
echo "  Visor interactivo:  $OUT/visor_mapa.html"
echo "  Genes en el mapa:   $OUT/genes_en_mapa.tsv"
echo "============================================================"
