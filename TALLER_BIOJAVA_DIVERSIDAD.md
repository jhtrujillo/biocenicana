# Taller de Bioinformática: Análisis de Diversidad Genética con BioJava

## Introducción
En el mejoramiento genético de cultivos poliploides como la caña de azúcar (*Saccharum spp.*), cuantificar la relación entre individuos es fundamental para la selección de parentales y el estudio de la estructura poblacional. Debido a la complejidad de los genomas poliploides, el cálculo de distancias genéticas debe considerar las dosis alélicas en lugar de genotipos discretos simples.

BioJava permite calcular matrices de distancia genética de forma eficiente utilizando modelos de Máxima Verosimilitud optimizados para poliploides.

---

## 1. Generación de la Matriz de Distancia Genética

Para este ejercicio, utilizaremos el set de datos filtrado de la población **CC 01-1940**.

### 0. Crear carpeta de resultados
Antes de empezar, crearemos una carpeta para organizar todos nuestros resultados:
```bash
mkdir -p taller_bioinformatica
```

### Opciones de Métodos de Distancia
BioJava permite elegir entre diferentes métricas matemáticas según el objetivo del estudio. Puedes probar cada una cambiando el parámetro `--method`:

```bash
# 1. Distancia de Nei (Ideal para filogenia y evolución)
java -jar target/biojava.jar genetic-distance -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method nei > taller_bioinformatica/matrix_nei.tsv

# 2. Distancia Euclidiana (Común en análisis de clusters)
java -jar target/biojava.jar genetic-distance -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method euclidean > taller_bioinformatica/matrix_euclidean.tsv

# 3. Distancia de Manhattan/p-distance (Default)
java -jar target/biojava.jar genetic-distance -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method manhattan > taller_bioinformatica/matrix_manhattan.tsv

# 4. Distancia de Rogers
java -jar target/biojava.jar genetic-distance -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method rogers > taller_bioinformatica/matrix_rogers.tsv

# 5. p-distance / IBS (Identity By State)
# Nota: En BioJava, p-distance e IBS son equivalentes a Manhattan (proporción de diferencias)
java -jar target/biojava.jar genetic-distance -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method ibs > taller_bioinformatica/matrix_ibs.tsv
```

### Explicación de parámetros:
- `-v`: Archivo VCF de entrada (datos de variantes).
- `-p`: Nivel de ploidía (10 para Saccharum).
- `--method`: Algoritmo matemático (`manhattan`, `euclidean`, `nei`, `rogers`, `ibs`, `p-distance`).
- `> taller_bioinformatica/archivo.tsv`: Redirección de la salida a la carpeta del taller.

> [!NOTE]
> **IBS y p-distance**: Miden la proporción de alelos compartidos o diferentes. Son métricas estándar para detectar duplicados o clones en una población.

---

## 2. Visualización Interactiva con Filogenia (Dashboard Web)

Una vez obtenida la matriz, el siguiente paso lógico es reconstruir un árbol filogenético para visualizar los clusters de forma interactiva en un navegador.

```bash
### Construcción del Árbol con diferentes Métricas
Al igual que con la matriz, puedes generar el árbol usando el método que prefieras. El visualizador web se adaptará a la métrica elegida:

```bash
# 1. Árbol basado en Distancia de Nei
java -jar target/biojava.jar snp-tree -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method nei -o taller_bioinformatica/arbol_nei.nwk

# 2. Árbol basado en Distancia Euclidiana
java -jar target/biojava.jar snp-tree -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method euclidean -o taller_bioinformatica/arbol_euclidean.nwk

# 3. Árbol basado en Manhattan (Default)
java -jar target/biojava.jar snp-tree -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method manhattan -o taller_bioinformatica/arbol_manhattan.nwk

# 4. Árbol basado en Distancia de Rogers
java -jar target/biojava.jar snp-tree -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method rogers -o taller_bioinformatica/arbol_rogers.nwk

# 5. Árbol basado en IBS / p-distance
java -jar target/biojava.jar snp-tree -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf -p 10 --method ibs -o taller_bioinformatica/arbol_ibs.nwk
```

*Cada ejecución generará un archivo `.html` correspondiente (ej: `arbol_nei.html`) que podrás abrir para explorar el árbol de forma interactiva.*

---

## 3. Análisis de Estructura Poblacional (PCA)

El Análisis de Componentes Principales (PCA) permite reducir la dimensionalidad de miles de SNPs a unos pocos ejes que explican la mayor parte de la variación genética. Esto es útil para detectar sub-poblaciones o grupos de parentesco.

### Comando de ejecución
Utiliza el subcomando `pop-structure`:

```bash
# Calcular PCA y Matriz de Parentesco (Kinship)
java -jar target/biojava.jar pop-structure \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -p 10 \
  -o taller_bioinformatica/pca_1940
```

### Resultados generados:
*   **`taller_bioinformatica/pca_1940.pca.csv`**: Tabla con las coordenadas de cada individuo en los primeros 10 componentes principales.
*   **`taller_bioinformatica/pca_1940.kinship.csv`**: Matriz de parentesco (VanRaden) usada para GWAS o selección genómica.
*   **`taller_bioinformatica/pca_1940.pca.html`**: **Visualizador web interactivo.** Permite rotar el gráfico en 3D y colorear por grupos detectados automáticamente.

---

---

## 4. Genómica Comparativa y Sintenia

La genómica comparativa nos permite entender cómo se organizan los genes entre diferentes genomas (sintenia) y cómo ha evolucionado la población en regiones genómicas específicas.

### Comando de ejecución
Utiliza el subcomando `comp-gen`. Este comando integra resultados de McScanX con datos de diversidad poblacional:

```bash
# 1. Crear carpeta para resultados si no existe
mkdir -p taller_bioinformatica

# 2. Integrar Sintenia, Diversidad Poblacional y Evolución (Ka/Ks)
java -jar target/biojava.jar comp-gen \
  --gff1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3 \
  --gff2 benchmarks/genomas/1940/CC-01-1940.gff3 \
  --collinearity benchmarks/genomica_comparativa/1940_vs_r570/mcscanx/1940_vs_r570.collinearity \
  --cds1 benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna \
  --cds2 benchmarks/genomas/1940/CC-01-1940.cds.fna \
  --vcf benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --viz taller_bioinformatica/visor_sintenia.html \
  --output taller_bioinformatica/reporte_comparativo.tsv \
  --name1 "R570" \
  --name2 "CC 1940" \
  --organism Saccharum
```

> [!TIP]
> **¿El gráfico sale vacío?** 
> Si al abrir `visor_sintenia.html` no ves nada, verifica la consola de Java. Si ves advertencias de `0 matched genes`, significa que los IDs en el archivo de colinealidad no coinciden con los del GFF. Asegúrate de estar usando los archivos GFF3 completos de la carpeta `benchmarks/genomas/`.

### ¿Qué estamos analizando?
1.  **Sintenia**: Identificación de bloques de genes que se mantienen en el mismo orden entre el genoma de referencia (R570) y el de estudio (CC 1940).
2.  **Densidad de SNPs**: El parámetro `--vcf` calcula cuántas variantes hay por cada bloque sinténico. Esto permite identificar regiones conservadas vs regiones hiper-variables.
3.  **Visualización**: El archivo `.html` genera un **Synteny Browser** interactivo donde puedes navegar por los cromosomas y ver los bloques de genes.
4.  **Análisis WGD (Ks Distribution)**: Al incluir `--cds1` y `--cds2`, BioJava calcula la tasa de sustitución sinónima (Ks). Esto permite identificar eventos de duplicación del genoma completo (*Whole Genome Duplication*) y estimar su antigüedad en millones de años.

### ⚠️ Resolución de Problemas (Sintenia)

Si ves el mensaje `[Warning] Block 0: 0 matched genes` o el visor HTML sale en blanco:

1.  **Nombres de Genes**: Asegúrate de que los IDs en el archivo `.collinearity` coincidan con el atributo `ID=` en el GFF3. BioJava ahora intenta limpiar prefijos como `gene:` o sufijos como `.1` automáticamente.
2.  **Orden de Genomas**: BioJava ahora detecta automáticamente si el archivo de colinealidad tiene los genomas invertidos respecto a los GFFs (e.g., G2 vs G1) y los intercambia en el reporte. Observa el log: `[Auto-Fix] Detected swapped genome order`.
3.  **Modo Offline**: El visor de sintenia (`visor_sintenia.html`) ahora incluye todas las librerías necesarias (D3.js) de forma interna. Debería funcionar directamente abriéndolo en Chrome o Safari sin necesidad de un servidor web.

Si el PCA no carga, asegúrate de tener conexión a internet ya que usa Plotly desde un CDN.

---

## 5. Reporte de Consenso de Parentesco

A menudo, los mejoradores necesitan cruzar la información de diferentes análisis para tomar una decisión final. El comando `rel-consensus` integra los resultados del PCA, la Distancia Genética y el Kinship en una sola tabla resumida.

### Comando de ejecución
```bash
java -jar target/biojava.jar rel-consensus \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -p 10 \
  -o taller_bioinformatica/reporte_consenso.csv
```

### ¿Qué obtenemos?
El comando genera un archivo `.csv` con las siguientes columnas para cada par de individuos:
1.  **Same_PCA_Cluster**: ¿Fueron agrupados juntos por el PCA?
2.  **Distance_PC_Space**: Distancia "física" en el mapa del PCA.
3.  **Kinship_VanRaden**: Valor de parentesco genómico real.
4.  **Inferred_Relationship**: Clasificación automática (Clon, Hermano, Medio Hermano, No relacionado).

### Interpretación Estratégica:
*   **Validación:** Si dos muestras están en el mismo clúster del PCA y su Kinship es > 0.40, puedes confirmar con total seguridad que son **familia directa**.
*   **Selección:** Para nuevos cruces, busca pares de muestras donde `Same_PCA_Cluster` sea **NO** y el Kinship sea cercano a **0** o negativo.

---

## 6. Interpretación de los Resultados

La salida es una **matriz N x N** (donde N es el número de muestras) separada por tabuladores.

### Ejemplo de Salida (Subset):
| | 100 | 101 | 102 |
|---|---|---|---|
| **100** | 0.0000 | 0.2451 | 0.3102 |
| **101** | 0.2451 | 0.0000 | 0.2894 |
| **102** | 0.3102 | 0.2894 | 0.0000 |

### Significado:
*   **Diagonal (0.0000)**: Indica la distancia de un individuo consigo mismo.
*   **Valores cercanos a 0**: Indican individuos genéticamente muy similares (posibles duplicados o clones).
*   **Valores altos**: Indican mayor divergencia genética.

---

## 5.5. Estimación de Ploidía por Individuo (Nuevo)

En poblaciones heterogéneas o cuando se trabaja con materiales de origen desconocido, la ploidía de cada individuo puede no ser uniforme. BioJava puede **estimar automáticamente** la ploidía de cada muestra directamente desde las frecuencias alélicas del VCF, sin necesidad de datos adicionales.

### ¿Cómo funciona?

El algoritmo realiza **dos pasadas** sobre el VCF:

1. **Pasada 1 — Estimación:** Para cada individuo y cada ploidía candidata `p`, calcula el error cuadrático medio entre las frecuencias alélicas observadas (BAF) y los niveles teóricos esperados (`0/p, 1/p, ... p/p`). La ploidía que produce el menor error es seleccionada.
2. **Pasada 2 — Dosaje:** Calcula la matriz de dosis alélicas usando la ploidía específica estimada para cada individuo.

> [!IMPORTANT]
> Si un individuo es **10x**, sus dosajes se calculan como `{0.0, 0.1, 0.2, ..., 1.0}`. Si otro es **8x**, sus dosajes serán `{0.0, 0.125, 0.25, ..., 1.0}`. Cada individuo usa su propio sistema de discretización.

### Comando de ejecución

```bash
# Estimar ploidía individualmente y calcular dosis alélicas
java -jar target/biojava.jar allele-dosage \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  --estimate-ploidy \
  --ploidy-candidates 2,4,6,8,10,12 \
  --ploidy-report taller_bioinformatica/reporte_ploidia.tsv \
  -i bsdp-mode \
  -md 5 \
  > taller_bioinformatica/dosis_por_ploidia.tsv
```

### Explicación de parámetros:
| Parámetro | Descripción |
|---|---|
| `--estimate-ploidy` | Activa la estimación automática. El parámetro `-p` se ignora como fuente principal. |
| `--ploidy-candidates` | Lista separada por comas de ploidías a evaluar. Acepta cualquier entero ≥ 2. |
| `--ploidy-report` | Archivo TSV donde se guarda el reporte de ploidías estimadas por individuo. |
| `-i bsdp-mode` | Usa conteos de lectura primero; imputa valores faltantes con la moda. |
| `-md 5` | Rechaza genotipos con menos de 5 lecturas totales (control de calidad). |
| `> dosis_por_ploidia.tsv` | Redirige la **matriz de dosaje** (stdout) al archivo de resultados. |

> [!NOTE]
> El reporte de ploidías se imprime en `stderr` (consola) si no se especifica `--ploidy-report`. La **matriz de dosaje** siempre va a `stdout`. Esto permite redirigirlos por separado.

### Ejemplo de reporte de ploidía generado (`reporte_ploidia.tsv`):

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║             [BioJava] Per-Individual Ploidy Estimation Report                   ║
╚══════════════════════════════════════════════════════════════════════════════════╝
Sample                              Ploidy_Estimated  Confidence  SNPs_Used    Status
─────────────────────────────────────────────────────────────────────────────────────
CC01-1940_replicate_1               10                0.921       18432        ✓ OK
CC01-1940_replicate_2               10                0.908       17891        ✓ OK
muestra_sospechosa_003              8                 0.073       12340        ⚠ AMBIGUOUS
─────────────────────────────────────────────────────────────────────────────────────
Total samples: 3  |  Ambiguous: 1
```

### Interpretación del reporte:
- **`Ploidy_Estimated`**: La ploidía que mejor explica el patrón de frecuencias alélicas del individuo.
- **`Confidence`** ∈ [0, 1]: Qué tan claro es el ajuste. `Confidence = 1 - (mejor_residual / segundo_mejor_residual)`. Valores altos indican señal clara.
- **`SNPs_Used`**: Número de SNPs con profundidad suficiente usados para la estimación. Valores bajos indican muestras con baja cobertura.
- **`⚠ AMBIGUOUS`**: Confianza < 0.10. El patrón BAF es igualmente compatible con dos ploidías. Revisa la cobertura o considera incluir esa ploidía como candidata adicional.

> [!TIP]
> Si obtienes muchas muestras `AMBIGUOUS`, intenta añadir ploidías impares a `--ploidy-candidates` (ej: `2,3,4,5,6,8,10,12`) o aumentar el filtro de profundidad con `-md 10`.

---

## 8. Anotación Funcional Avanzada y Matriz de Haplotipos

Esta es la herramienta más completa de BioJava. Permite integrar en un solo tablero interactivo la función de los genes, el efecto de las mutaciones (VEP), la recomendación de marcadores para laboratorio (KASP) y la estructura de haplotipos de la población.

### Comando de ejecución
Utiliza el subcomando `annotate`. Este comando requiere el genoma de referencia para el análisis KASP:

```bash
# Generar la Suite Profesional de Anotación y Población
java -jar target/biojava.jar annotate \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -g benchmarks/genomas/1940/CC-01-1940.gff3 \
  -r benchmarks/genomas/1940/CC-01-1940.fasta \
  --gff2 benchmarks/genomas/R570/R570.gff3 \
  -p benchmarks/genomas/1940/CC-01-1940.protein.faa \
  -c benchmarks/genomas/1940/CC-01-1940.cds.fna \
  -w 5000 \
  -o taller_bioinformatica/BioJava_Suite_Final.html
```

### ¿Qué aprenderemos en este módulo?

1.  **Buscador Híbrido**: Usa la barra lateral para buscar por funciones (ej. "sugar", "stress", "transporter"). El sistema autocompletará términos de GO e InterPro.
2.  **Variant Effect Predictor (VEP)**: Identifica si un SNP cambia un aminoácido (**Missense**) o si es silencioso (**Synonymous**).
3.  **Recomendación KASP (⭐⭐⭐)**: 
    *   Busca los SNPs con **3 estrellas**. Son los mejores candidatos para el laboratorio porque no tienen interferencia de otros SNPs cercanos y tienen un GC ideal.
    *   Haz clic en **"KASP Primer"** para obtener la secuencia lista para pedir a tu proveedor.
4.  **Matriz de Haplotipos**: 
    *   Haz clic en el botón **"Haplotypes"**.
    *   Observa cómo se distribuyen los genotipos (Verde/Naranja/Rojo) en tus 220 muestras para todos los genes que te interesan. ¿Ves bloques de herencia?

---

## 9. Catálogo Completo de Comandos (BioJava Toolkit)

BioJava es una navaja suiza para genómica de poliploides. Aquí tienes todos los comandos disponibles que puedes usar en tu taller:

### A. Gestión y Control de Calidad (QC)
*   **`vcf-stats`**: Genera reportes estadísticos y tableros interactivos de calidad.
*   **`vcf-filter`**: Filtra variantes por MAF, missingness y equilibrio de Hardy-Weinberg.
*   **`vcf-merge`**: Combina múltiples archivos VCF en uno solo.

### B. Análisis de Diversidad y Estructura
*   **`pop-structure`**: Realiza Análisis de Componentes Principales (PCA) para ver agrupamientos.
*   **`genetic-distance`**: Calcula la matriz de distancia genética.
*   **`snp-tree`**: Reconstruye árboles filogenéticos Neighbor-Joining.
*   **`snp-explorer`**: Permite auditar visualmente el comportamiento de cada SNP.
*   **`rel-consensus`**: Genera el reporte integrado de parentesco y grupos.

### C. Anotación y Genómica Funcional (Suite Profesional)
*   **`annotate`**: **[NUEVO]** Genera el tablero integral con VEP, KASP, Sintenia y Matriz de Haplotipos.
*   **`allele-dosage`**: Exporta las dosis alélicas puras (0, 1, 2... k).
*   **`kaks-calc`**: Calcula tasas de sustitución Ka/Ks para detectar selección.
*   **`comp-gen`**: Integra sintenia, WGD y anotaciones funcionales globales.

### D. Exportación a otros Formatos
*   **`gwaspoly-export`**: Prepara datos para el paquete R GWASpoly.
*   **`joinmap`**: Convierte datos para mapeo de ligamiento en JoinMap.

---

## Conclusión
El uso de la Suite BioJava permite a los mejoradores capturar la verdadera variación genética en poliploides, facilitando decisiones más precisas desde la estructura poblacional hasta la validación de marcadores funcionales en el laboratorio.
