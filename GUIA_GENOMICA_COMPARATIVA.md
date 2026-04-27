# Guía Completa: Genómica Comparativa con BioCenicana
## Caso de estudio: *Saccharum hybrid* R570 vs CC 01-1940

Esta guía lleva al usuario desde cero hasta un análisis genómico comparativo multicapa,
usando datos reales de caña de azúcar disponibles en el directorio `benchmarks/`.

---

## 0. Preparación: Compilar y Empaquetar el JAR

Este paso solo se necesita hacer **una vez** (o cada vez que el código cambia).

```bash
mvn package -q
```

A partir de aquí, todos los comandos usan:

```bash
java -jar target/biocenicana.jar <subcomando> [opciones]
```

Puedes verificar que todo funciona con:

```bash
java -jar target/biocenicana.jar --help
```

---

## 1. Archivos de Entrada

| Componente | Ruta |
| :--- | :--- |
| **GFF3 — R570** | `benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3` |
| **GFF3 — 1940** | `benchmarks/genomas/1940/CC-01-1940.gff3` |
| **Colinealidad (McScanX)** | `benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity` |
| **CDS — R570** | `benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna` |
| **CDS — 1940** | `benchmarks/genomas/1940/CC-01-1940.cds.fna` |
| **VCF — 1940** | `benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf` |

---

## Paso 1 — Sintenia y Diversidad Genómica

Identifica los bloques sinténicos entre R570 y 1940, y calcula la densidad de SNPs por bloque usando el VCF real.

```bash
java -jar target/biocenicana.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --viz results/paso1_sintenia_diversidad.html
```

**Resultados:**
- `results/paso1_sintenia_diversidad.html` — Dashboard interactivo con mapa de calor de diversidad.
- `comp_gen_report.tsv` — Reporte tabular con todos los pares de genes sinténicos.

---

## Paso 2 — Anotación Funcional y Enriquecimiento GO

Conecta cada bloque sinténico con sus funciones biológicas. Identifica qué procesos biológicos están enriquecidos en cada bloque (Test Exacto de Fisher, p < 0.05).

```bash
java -jar target/biocenicana.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --annot1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --annot2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --viz results/paso2_anotaciones_go.html
```

**Resultados:**
- `results/paso2_anotaciones_go.html` — Al pasar el cursor sobre cada bloque, el tooltip muestra la función de cada gen y los términos GO enriquecidos.

---

## Paso 3 — Detección de WGD y Exportación de Ortólogos para Filogenia

Integra las secuencias CDS para:
1. Calcular **Ks** (sustituciones sinónimas) y detectar eventos de duplicación del genoma (WGD).
2. Identificar pares de ortólogos estrictos 1:1, alinearlos (Needleman-Wunsch) y exportar la **super-matriz filogenética**.

```bash
java -jar target/biocenicana.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --cds1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --annot1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --annot2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --export-orthologs results/supermatrix_1940_r570.fasta \
  --viz results/paso3_wgd_filogenia.html
```

**Resultados:**
- `results/paso3_wgd_filogenia.html` — Dashboard con histograma de distribución Ks integrado.
- `results/supermatrix_1940_r570.fasta` — Super-matriz FASTA lista para usar en IQ-TREE o RAxML.

---

## Paso 4 — Dashboard Final Integrado (Todas las Capas)

Consolida todos los análisis anteriores en un único explorador interactivo.

```bash
java -jar target/biocenicana.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --cds1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --annot1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --annot2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --vcf benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --export-orthologs results/supermatrix_1940_r570.fasta \
  --viz results/dashboard_final.html
```

**Resultados finales:**
- `results/dashboard_final.html` — Explorador completo: Ribbons, Dotplot 2D, Vista Circos, filtros de diversidad/Ka/Ks, enriquecimiento GO y detección de SVs.
- `results/supermatrix_1940_r570.fasta` — Insumo para análisis filogenético.
- `comp_gen_report.tsv` — Reporte TSV con todas las métricas integradas.

---

## Guía del Dashboard Interactivo

| Modo de Vista | Uso Principal |
| :--- | :--- |
| **Ribbons** | Ver la correspondencia directa entre cromosomas |
| **Dotplot 2D** | Identificar inversiones y translocaciones a gran escala |
| **Circos** | Resumen global de toda la comparación genómica |

| Filtro | Función |
| :--- | :--- |
| **Block Size** | Elimina bloques pequeños (ruido) |
| **Diversity Mode** | Colorea bloques por densidad de SNPs |
| **Ka/Ks Mode** | Muestra presión de selección: azul = purificadora, rojo = positiva |
| **Búsqueda de Gen** | Localiza cualquier gen y ve su contexto sinténico |
