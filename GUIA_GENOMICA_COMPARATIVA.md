# Guía Maestra: Análisis de Genómica Comparativa con BioCenicana

Esta guía utiliza datos reales del proyecto BioCenicana (Cultivar R570 vs CC 01-1940) para demostrar el flujo de trabajo completo del módulo `comp-gen`.

## 1. Descripción del Dataset de Prueba
Para este tutorial utilizaremos los archivos ubicados en `benchmarks/`:

| Componente | Ruta del Archivo |
| :--- | :--- |
| **GFF3 (R570)** | `benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3` |
| **GFF3 (1940)** | `benchmarks/genomas/1940/CC-01-1940.gff3` |
| **Colinealidad** | `benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity` |
| **CDS (R570)** | `benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna` |
| **CDS (1940)** | `benchmarks/genomas/1940/CC-01-1940.cds.fna` |
| **Variantes (VCF)** | `benchmarks/vcfs/maize/maize.vcf` (Para pruebas de densidad de SNPs) |

---

## 2. Ejecución del Análisis Completo

El siguiente comando integra todas las capas analíticas: sintenia, filogenia, variantes estructurales y funciones GO.

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
  --export-orthologs results/supermatrix_1940_r570.fasta \
  --viz results/dashboard_interactivo.html"
```

---

## 3. Funciones Analíticas Detalladas

### A. Detección de Eventos de Duplicación (WGD)
Al incluir los archivos CDS, BioCenicana calcula automáticamente los valores de **Ks (Sustituciones Sinónimas)** para cada par de genes.
- **Resultado:** En el dashboard interactivo, verás un histograma de frecuencias de Ks.
- **Uso:** Los picos en valores bajos de Ks indican duplicaciones recientes, mientras que picos en valores altos indican eventos ancestrales de poliploidización.

### B. Exportación de Super-Matriz de Ortólogos
La opción `--export-orthologs` identifica genes que mantienen una relación 1:1 estricta dentro de los bloques sinténicos.
- **Proceso:** El sistema alinea cada par usando el algoritmo Needleman-Wunsch y concatena los alineamientos.
- **Resultado:** Un archivo FASTA (`supermatrix_1940_r570.fasta`) listo para inferencia filogenética.

### C. Integración de Variantes Estructurales (SVs)
Si posees un archivo VCF con variantes estructurales (obtenido mediante herramientas como Delly o Lumpy), puedes usar `--sv`.
- **Detección:** El sistema busca solapamientos espaciales entre los bloques de sintenia y las coordenadas de los SVs.
- **Visualización:** Los bloques afectados se marcan con una advertencia `⚠ SV Detected` y un estilo visual diferenciado (bordes punteados).

### D. Enriquecimiento Funcional (GO)
Usando los atributos del GFF3 o un archivo de anotación separado, BioCenicana realiza un **Test Exacto de Fisher** por cada bloque.
- **Propósito:** Identificar si un bloque sinténico contiene una sobre-representación de genes relacionados con una función biológica específica (ej. "respuesta a estrés hídrico").
- **Visualización:** Aparece en el tooltip al pasar el cursor sobre cualquier bloque en el explorador.

---

## 4. Guía de Navegación del Dashboard

El archivo `dashboard_interactivo.html` permite tres modos de visualización:

1.  **Vista de Cintas (Ribbons):** Ideal para ver la correspondencia directa entre cromosomas. El color de las cintas puede representar el ID del bloque o el valor de Ka/Ks (presión de selección).
2.  **Dotplot 2D:** Esencial para identificar inversiones a gran escala y translocaciones complejas entre los dos cultivares de caña.
3.  **Vista Circos:** Proporciona un resumen global de todas las relaciones sinténicas en un formato circular compacto.

### Filtros en el Panel Lateral:
- **Block Size:** Filtra bloques pequeños para limpiar el ruido.
- **Orientation:** Separa bloques colineales de bloques invertidos.
- **Diversity Mode:** Cambia el color de los bloques para mostrar la densidad de SNPs detectada en el VCF.

---

## 5. Solución de Problemas Comunes

- **Error de Memoria:** Si los archivos BLAST/Collinearity son muy grandes, aumenta la memoria de Maven con `MAVEN_OPTS="-Xmx8G"`.
- **Cromosomas no encontrados:** Asegúrate de que los nombres de los cromosomas en el archivo de colinealidad coincidan exactamente con los del GFF.
- **Alineamiento Lento:** La exportación de ortólogos puede tardar varios minutos si hay miles de pares 1:1; el sistema mostrará el progreso en la consola.
