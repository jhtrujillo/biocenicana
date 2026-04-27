# Tutorial: Genómica Comparativa Avanzada con BioCenicana

Este tutorial describe cómo utilizar el módulo `comp-gen` de BioCenicana para integrar resultados de sintenia (McScanX/SynMap), datos de variantes estructurales, selección evolutiva y anotaciones funcionales.

## 1. Requisitos Previos

Para este análisis, utilizaremos los datos de referencia ubicados en el directorio `benchmarks/`:
- **GFF3**: Archivos de anotación de los genomas comparados.
- **Collinearity**: Salida de McScanX que define los bloques sinténicos.
- **CDS FASTA**: Secuencias codificantes para alineamiento de ortólogos.
- **VCF**: Datos de variantes para diversidad poblacional y variantes estructurales (SV).

## 2. Estructura del Comando

El comando principal es `comp-gen`. Aquí se detallan los parámetros clave:

### Parámetros de Estructura (Obligatorios)
- `--gff1` / `--gff2`: Archivos GFF3 de los genomas.
- `--collinearity`: Archivo `.collinearity` generado por McScanX.

### Parámetros de Evolución y Función (Nuevos)
- `--cds1` / `--cds2`: Secuencias CDS para extraer longitudes y realizar alineamientos 1:1.
- `--annot1` / `--annot2`: Archivos funcionales (pueden ser los mismos GFF3).
- `--vcf`: Archivo VCF para calcular densidad de SNPs por bloque.
- `--sv`: Archivo VCF con variantes estructurales (DELETION, INVERSION).
- `--kaks`: Archivo de resultados de Ka/Ks para visualizar presión de selección.
- `--export-orthologs`: Ruta para exportar la super-matriz de ortólogos alineados.
- `--viz`: Ruta del dashboard interactivo HTML.

## 3. Ejemplo Práctico: Caña de Azúcar (R570 vs CC 01-1940)

Ejecuta el siguiente comando para realizar un análisis completo integrando todas las capas analíticas:

```bash
mvn exec:java -Dexec.mainClass="org.cenicana.bio.Main" -Dexec.args="comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --cds1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --annot1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --annot2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --vcf benchmarks/vcfs/maize/maize.vcf \
  --export-orthologs phylogeny_supermatrix.fasta \
  --viz synteny_explorer.html"
```

## 4. Interpretación de Resultados

### A. Reporte Integrado (TSV)
El archivo `comp_gen_report.tsv` contiene una fila por cada par de genes, incluyendo:
- IDs de genes y coordenadas.
- E-value del alineamiento original.
- Longitudes de CDS y proteínas.
- Estado (Syntenic vs Orphan).

### B. Super-matriz de Ortólogos (FASTA)
El archivo `phylogeny_supermatrix.fasta` contiene las secuencias de todos los ortólogos 1:1 alineados y concatenados. Este archivo está listo para ser usado en herramientas como `RAxML`, `IQ-TREE` o el módulo `phylogeny` de BioCenicana.

### C. Explorador Interactivo (HTML)
Abre `synteny_explorer.html` en tu navegador para acceder a las 4 capas de análisis:

1.  **Dashboard de Duplicación (WGD):** Observa el histograma de Ks para identificar picos evolutivos que corresponden a eventos de duplicación del genoma completo.
2.  **Mapa de Calor Evolutivo:** Los bloques cambian de color según la densidad de SNPs o los valores de Ka/Ks (Selección Purificadora vs Positiva).
3.  **Detección de SVs:** Los bloques interrumpidos por inversiones o grandes deleciones se resaltan con bordes punteados.
4.  **Enriquecimiento GO:** Al pasar el cursor sobre un bloque, verás los términos de Gene Ontology significativamente enriquecidos (p < 0.05).

## 5. Tips Avanzados

- **Filtros Dinámicos:** En el sidebar del dashboard, puedes filtrar bloques por número mínimo de genes para eliminar "ruido" de alineamientos cortos.
- **Búsqueda de Genes:** Usa el buscador para localizar rápidamente un gen de interés y ver su contexto sinténico en el otro genoma.
- **Modos de Vista:** Cambia entre **Ribbons** (vista clásica), **Dotplot 2D** (para detectar inversiones masivas) y **Circos** (para una visión genómica global).
