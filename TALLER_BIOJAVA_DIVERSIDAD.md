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

## 7. Catálogo Completo de Comandos (BioJava Toolkit)

BioJava es una navaja suiza para genómica de poliploides. Aquí tienes todos los comandos disponibles que puedes usar en tu taller:

### A. Gestión y Control de Calidad (QC)
*   **`vcf-stats`**: Genera reportes estadísticos y tableros interactivos de calidad.
    ```bash
    java -jar target/biojava.jar vcf-stats -v input.vcf -o reporte_qc
    ```
*   **`vcf-filter`**: Filtra variantes por MAF, missingness y equilibrio de Hardy-Weinberg.
    ```bash
    java -jar target/biojava.jar vcf-filter -v input.vcf --min-maf 0.05 --max-missing 0.2 -o filtrado.vcf
    ```
*   **`vcf-merge`**: Combina múltiples archivos VCF en uno solo.
    ```bash
    java -jar target/biojava.jar vcf-merge -i batch1.vcf,batch2.vcf -o unificado.vcf
    ```

### B. Análisis de Diversidad y Estructura
*   **`pop-structure`**: Realiza Análisis de Componentes Principales (PCA) para ver agrupamientos.
    ```bash
    java -jar target/biojava.jar pop-structure -v filtrado.vcf -p 10 -o estructura
    ```
*   **`genetic-distance`**: Calcula la matriz de distancia genética.
*   **`snp-tree`**: Reconstruye árboles filogenéticos Neighbor-Joining.
    ```bash
    java -jar target/biojava.jar snp-tree -v filtrado.vcf -p 10 -o arbol.nwk
    ```
*   **`snp-explorer`**: Permite auditar visualmente el comportamiento de cada SNP.
    ```bash
    java -jar target/biojava.jar snp-explorer --vcf filtrado.vcf --pca estructura.pca.csv -o visor_snps.html
    ```
*   **`rel-consensus`**: Genera el reporte integrado de parentesco y grupos.
    ```bash
    java -jar target/biojava.jar rel-consensus -v filtrado.vcf -p 10 -o reporte.csv
    ```

### C. Desequilibrio de Ligamiento y Genómica Funcional
*   **`ld`**: Calcula el decaimiento de LD (Linkage Disequilibrium).
    ```bash
    java -jar target/biojava.jar ld -v filtrado.vcf -o reporte_ld
    ```
*   **`allele-dosage`**: Exporta las dosis alélicas puras (0, 1, 2... k).
    ```bash
    java -jar target/biojava.jar allele-dosage -v filtrado.vcf -p 10 > dosis.tsv
    ```

### D. Genómica Comparativa y Evolución
*   **`comp-gen`**: Integra sintenia, WGD y anotaciones funcionales.
    ```bash
    java -jar target/biojava.jar comp-gen --gff1 g1.gff --gff2 g2.gff --collinearity col.txt --viz sintenia.html
    ```
*   **`kaks-calc`**: Calcula tasas de sustitución Ka/Ks para detectar selección.
    ```bash
    java -jar target/biojava.jar kaks-calc --cds1 g1.fasta --cds2 g2.fasta -o seleccion.tsv
    ```

### E. Exportación a otros Formatos
*   **`gwaspoly-export`**: Prepara datos para el paquete R GWASpoly.
*   **`joinmap`**: Convierte datos para mapeo de ligamiento en JoinMap.

---

## Conclusión
El uso de matrices de distancia basadas en dosis alélicas permite a los mejoradores capturar la verdadera variación genética en poliploides, facilitando decisiones más precisas en los programas de cruzamiento.
