#!/usr/bin/env bash
# =============================================================================
# Pipeline: Construcción del mapa genético y ubicación de genes de sacarosa
# Actividad 4.1.5 — Asociación de genes a marcadores/regiones de interés
#
# Fases:
#   1. Compilar BioJava
#   2. Construir mapa genético (BioJava genetic-map)
#   3. Posicionar genes de sacarosa en el mapa (map_genes_to_map.py)
#   4. Análisis de LD y tag SNPs por gen candidato (ld_candidate_genes.py)
#
# Uso:
#   bash run_genetic_map_pipeline.sh
#
# Requisitos:
#   - Java 11+ con Maven compilado (target/biojava.jar)
#   - Python 3 con entorno virtual (.venv)
#   - VCF de la población biparental en benchmarks/vcfs/mapa_genetico/
#   - VCF de variantes CC 01-1940 en benchmarks/vcfs/1940/
#   - GFF3 de CC 01-1940 en benchmarks/genomas/1940/
# =============================================================================

set -e

echo "============================================================"
echo "  PIPELINE: MAPA GENÉTICO + GENES DE SACAROSA              "
echo "============================================================"

# ── Rutas ────────────────────────────────────────────────────────────────────
VCF_MAP="benchmarks/vcfs/mapa_genetico/AllSamples_variants_geneticmap_pseudochromosomes_standarfilters_minInd_85_single_dosage.vcf"
VCF_ANN="benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf"
GFF="benchmarks/genomas/1940/CC-01-1940.gff3"
OUT="genomica_comparativa/mapa_genetico"
PYTHON=".venv/bin/python"

# Parámetros del mapa
PLOIDY=10
LOD=3.0
MAX_R=0.35
THIN_KB=500
MIN_LG=3
MAPPING_FUNC="kosambi"

# Parámetros LD (Fase 4)
LD_WINDOW_KB=500   # ventana a cada lado del gen para buscar marcadores

mkdir -p "$OUT"

# ── Paso 1: Compilar BioJava ─────────────────────────────────────────────────
echo ""
echo "[Paso 1/3] Compilando BioJava..."
mvn clean package -DskipTests -q
echo "  ✓ target/biojava.jar generado"

# ── Paso 2: Construir mapa genético ─────────────────────────────────────────
echo ""
echo "[Paso 2/3] Construyendo mapa genético..."
echo "  VCF:          $VCF"
echo "  Ploidía:      $PLOIDY"
echo "  LOD mínimo:   $LOD"
echo "  Max recomb:   $MAX_R"
echo "  Raleo:        $THIN_KB kb"
echo "  Min marcadores/LG: $MIN_LG"
echo ""

java -Xmx8g -jar target/biojava.jar genetic-map \
  -i "$VCF_MAP" \
  -p $PLOIDY \
  -o "$OUT/mapa_biparental_final.map" \
  --pseudo-only \
  --thin-kb $THIN_KB \
  --min-lg-markers $MIN_LG \
  --lod $LOD \
  --max-r $MAX_R \
  --mapping-function $MAPPING_FUNC \
  --viz "$OUT/visor_mapa_final.html"

echo ""
echo "  ✓ Mapa guardado en: $OUT/mapa_biparental_final.map"
echo "  ✓ Visor HTML:       $OUT/visor_mapa_final.html"

# ── Paso 3: Posicionar genes de sacarosa en el mapa ─────────────────────────
echo ""
echo "[Paso 3/3] Posicionando genes de sacarosa en el mapa..."

$PYTHON scripts/map_genes_to_map.py \
  --map "$OUT/mapa_biparental_final.map" \
  --gff "$GFF" \
  --output "$OUT/genes_en_mapa.tsv"

echo ""
echo "  ✓ Genes posicionados en: $OUT/genes_en_mapa.tsv"

# ── Paso 4: LD y tag SNPs por gen candidato ──────────────────────────────────
echo ""
echo "[Paso 4/4] Calculando LD y tag SNPs por gen candidato..."
echo "  VCF variantes: $VCF_ANN"
echo "  Ventana:       ±${LD_WINDOW_KB} kb"
echo ""

$PYTHON scripts/ld_candidate_genes.py \
  --vcf "$VCF_ANN" \
  --out-dir "$OUT/fase4_ld" \
  --window-kb $LD_WINDOW_KB

echo ""
echo "  ✓ Resumen LD:  $OUT/fase4_ld/ld_summary.tsv"
echo "  ✓ Tag SNPs:    $OUT/fase4_ld/tag_snps.tsv"

# ── Resumen final ────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  PIPELINE COMPLETADO — Actividad 4.1.5                    "
echo "============================================================"
echo ""
echo "  Fase 1 — Compilación:"
echo "    target/biojava.jar"
echo ""
echo "  Fase 2 — Mapa genético:"
echo "    $OUT/mapa_biparental_final.map"
echo "    $OUT/visor_mapa_final.html   (abrir en navegador)"
echo ""
echo "  Fase 3 — Genes de sacarosa en el mapa:"
echo "    $OUT/genes_en_mapa.tsv"
echo ""
echo "  Fase 4 — LD y marcadores por gen candidato:"
echo "    $OUT/fase4_ld/ld_summary.tsv"
echo "    $OUT/fase4_ld/tag_snps.tsv"
echo ""
echo "  Para ajustar parámetros edita las variables al inicio de este script."
echo ""
