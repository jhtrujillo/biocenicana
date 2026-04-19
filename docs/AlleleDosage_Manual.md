# 🧬 Manual de Usuario: Allele Dosage Calculator (Biocenicana)

Este manual te guiará paso a paso en el uso de la herramienta `AlleleDosageCalculator`, comenzando desde el uso más básico y tradicional, hasta llegar a las opciones avanzadas impulsadas por Inteligencia Artificial y Estadística. En cada nivel encontrarás ejemplos prácticos para entender exactamente qué le pasa a tus datos.

---

## 📖 1. El Contexto Biológico: ¿Qué es la "Dosis Alélica"?

En organismos **diploides** (como los humanos), heredamos 2 copias de cada gen. Por lo tanto, si miramos una mutación (SNP), solo hay 3 posibilidades:
- `0.0`: Tienes 0 copias de la mutación.
- `0.5`: Tienes 1 copia de la mutación (Heterocigoto).
- `1.0`: Tienes 2 copias de la mutación.

Sin embargo, organismos como la **caña de azúcar, la papa o el trigo** son **poliploides**. La caña de azúcar, por ejemplo, ¡puede tener hasta 10 copias de cada cromosoma! 
En un organismo con ploidía 10, un individuo puede tener 0, 1, 2... o 10 copias de esa mutación.

A la proporción de copias mutadas se le llama **Dosis Alélica**. La calculamos dividiendo las lecturas que capturó el secuenciador:
`Dosis = Lecturas de Referencia / (Lecturas de Referencia + Lecturas Alternas)`

---

## 🟢 Nivel 1: El Modo Clásico (Básico)

Si solo quieres calcular las dosis matemáticas puras y duras (el comportamiento clásico del programa), no necesitas usar ninguna bandera avanzada. El programa leerá tu archivo VCF, calculará la división matemática, forzará el resultado a encajar en el nivel de ploidía más cercano.

**💻 Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 --impute mode
```

**🔬 Ejemplo VCF en acción (Ploidía = 4):**
Imagina que el VCF tiene el formato GATK `GT:AD` donde `AD` es `LecturasRef,LecturasAlt`.
> **VCF Entrada:** `Individuo_1` tiene `0/1:20,80` (20 Ref, 80 Alt).
> **Proceso:** Dosis cruda = 20 / (20 + 80) = 0.20.
> **Salida:** En un tetraploide, los niveles son (0, 0.25, 0.50, 0.75, 1.0). El 0.20 se redondea a **0.25**.

---

## 🟡 Nivel 2: Filtro Estadístico de Calidad (`--min-depth`)

### El Problema (Ruido de Secuenciación)
Si el secuenciador solo capturó **2 fragmentos de ADN**, la matemática dirá: `1 / (1 + 1) = 0.5`. ¡Pero la planta tiene 10 copias! Confiar en ese dato es **estadísticamente peligroso**.

### La Solución
El parámetro `--min-depth` (o `-md`) te permite exigir un mínimo de "profundidad" de lecturas. 

**💻 Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 --impute mode -md 15
```

**🔬 Ejemplo VCF en acción:**
> **VCF Entrada:** 
> `Individuo_A`: `0/1:2,1` (Total lecturas: 3)
> `Individuo_B`: `0/1:40,40` (Total lecturas: 80)
> **Proceso:** El programa evalúa el `-md 15`. 
> **Salida:** `Individuo_B` supera el filtro y se le asigna dosis **0.50**. `Individuo_A` es descartado por falta de calidad y se marca como **"Dato Faltante"** (y luego se imputará con la moda).

---

## 🟠 Nivel 3: Imputación por Inteligencia Artificial (`--impute knn`)

### El Problema (Imputación Clásica Deficiente)
Si usas `--impute mode` para rellenar el dato faltante del `Individuo_A` (del ejemplo anterior), estás copiando lo que haga la mayoría de la población. Si la población mezcla familias distintas, arruinarás la matriz genética.

### La Solución (K-Nearest Neighbors)
Activando el método `knn`, el programa analiza el genoma entero y construye una matriz de parentesco genético para imputar usando solo a la familia de la muestra.

**💻 Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 -md 15 --impute bsdp-knn -k 2
```

**🔬 Ejemplo VCF en acción:**
> **Contexto:** Al `Individuo_A` le falta el dato del SNP_1. El programa busca a sus **2 vecinos genéticamente más parecidos** (`-k 2`) en todo el archivo VCF. Descubre que sus "hermanos" son el `Individuo_H1` y el `Individuo_H2`.
> **Proceso:** `Individuo_H1` tiene dosis 0.25 en el SNP_1. `Individuo_H2` tiene dosis 0.25 en el SNP_1.
> **Salida:** El programa rellena el hueco del `Individuo_A` con **0.25**, respetando su herencia genética, sin importar lo que opine el resto de la población.

---

## 🔴 Nivel 4: Redondeo Adaptativo por Machine Learning (`--adaptive-rounding`)

### El Problema (Reference Bias)
Los secuenciadores reales leen mejor el ADN original que el mutado. Esto causa que un individuo que debería tener una dosis de `0.25` termine reportando `0.38`. El redondeo matemático clásico se equivocaría y le asignaría `0.50`.

### La Solución (Clustering K-Means 1D)
El parámetro `--adaptive-rounding` (o `-ar`) enciende un algoritmo de Machine Learning No Supervisado que mira la nube poblacional para encontrar los verdaderos grupos (clústeres).

**💻 Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 --adaptive-rounding
```

**🔬 Ejemplo VCF en acción:**
> **VCF Entrada (SNP_1):** 
> `Individuo_1`: `0/1:36,64` (Frecuencia cruda = 0.36)
> `Individuo_2`: `0/1:38,62` (Frecuencia cruda = 0.38)
> `Individuo_3`: `0/1:37,63` (Frecuencia cruda = 0.37)
> **Proceso (Clásico vs IA):** 
> * Modo Clásico: Redondearía matemáticamente a **0.50** (¡Error!).
> * Modo IA (`-ar`): El clustering nota que *todo el grupo poblacional* está desplazado. Deduce que `0.37` es en realidad el nuevo centro desviado del grupo biológico `0.25`.
> **Salida:** Les asigna a los tres individuos la dosis corregida de **0.25**.

---

## 🚀 Nivel 5: El Modo "Dios" (Todo en Uno)

Al combinar todas estas funciones, obtienes el cálculo de dosis más robusto y estadísticamente avanzado posible en bioinformática de poliploides:

**💻 Comando Maestro:**
```bash
java -jar biocenicana.jar allele-dosage \
    -v mis_datos_crudos.vcf \
    -p 10 \
    -md 15 \
    --adaptive-rounding \
    --impute bsdp-knn \
    -k 7
```

**Resumen de la magia biológica que ocurre aquí:**
1. Elimina cualquier dato basura con menos de 15 lecturas (`-md 15`).
2. Usa clustering poblacional para esquivar los sesgos del secuenciador y asignar dosis exactas (`--adaptive-rounding`).
3. A los huecos dejados en el paso 1, los reconstruye copiando el ADN de sus 7 familiares más parecidos usando Inteligencia Artificial KNN (`--impute bsdp-knn -k 7`).
