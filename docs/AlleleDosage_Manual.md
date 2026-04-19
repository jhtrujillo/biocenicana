# 🧬 Manual de Usuario: Allele Dosage Calculator (Biocenicana)

Este manual te guiará paso a paso en el uso de la herramienta `AlleleDosageCalculator`, comenzando desde el uso más básico y tradicional, hasta llegar a las opciones avanzadas impulsadas por Inteligencia Artificial y Estadística.

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

Si solo quieres calcular las dosis matemáticas puras y duras (el comportamiento clásico del programa), no necesitas usar ninguna bandera avanzada. 

El programa leerá tu archivo VCF, calculará la división matemática, forzará el resultado a encajar en el nivel de ploidía más cercano, y rellenará los datos vacíos con la "moda" (el valor más común en la población).

**💻 Ejemplo de Comando Básico:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 10 --impute mode
```
* **Explicación bioinformática:** Calcula la dosis para una ploidía de 10 (`-p 10`). Si a una muestra le faltan datos (`./.`), rellénalo con el valor que más se repite en el resto de individuos (`--impute mode`).

---

## 🟡 Nivel 2: Filtro Estadístico de Calidad (`--min-depth`)

### El Problema (Ruido de Secuenciación)
Si el secuenciador solo capturó **2 fragmentos de ADN** para un individuo (uno de Referencia y uno Alterno), la matemática dirá: `1 / (1 + 1) = 0.5`. ¡Pero la planta tiene 10 copias! Es imposible estar seguros de que su dosis real es `0.5` con solo 2 lecturas. Confiar en ese dato es **estadísticamente peligroso**.

### La Solución
El parámetro `--min-depth` (o `-md`) te permite exigir un mínimo de "profundidad" de lecturas. Si un individuo no alcanza ese mínimo, el programa descarta sus lecturas y lo marca como un "dato faltante" (para imputarlo después).

**💻 Ejemplo de Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 10 --impute mode -md 15
```
* **Explicación bioinformática:** Exige que cada individuo tenga al menos 15 lecturas de ADN en ese SNP (`-md 15`). Si tiene 14 o menos, se marca como vacío y se rellena con la moda (`--impute mode`).

---

## 🟠 Nivel 3: Imputación por Inteligencia Artificial (`--impute knn`)

### El Problema (Imputación Clásica Deficiente)
Si usas `--impute mode` o `--impute mean` para rellenar los datos faltantes (o los que eliminaste en el Nivel 2), estás mezclando a toda la población. Si tienes plantas de dos familias distintas, promediarlas destruirá la estructura genética de tus datos.

### La Solución (K-Nearest Neighbors)
Activando el método `knn`, el programa analiza el genoma entero y construye una matriz de parentesco genético. Cuando encuentra un dato faltante, busca a los **K individuos genéticamente más parecidos** (sus primos/vecinos) y usa la dosis de ellos para rellenar el hueco.

**💻 Ejemplo de Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 10 -md 15 --impute bsdp-knn -k 5
```
* **Explicación bioinformática:** Imputa usando KNN (`--impute bsdp-knn`). Busca a los 5 familiares más cercanos (`-k 5`) de la muestra que tiene el dato faltante para promediar el valor.

---

## 🔴 Nivel 4: Redondeo Adaptativo por Machine Learning (`--adaptive-rounding`)

### El Problema (Reference Bias)
Los secuenciadores leen mejor el ADN original que el mutado. Esto causa que un individuo que debería ser `0.25` termine reportando `0.38`. El método matemático estricto se confundiría y diría: *"0.38 está más cerca de 0.50, le pondré 0.50"*. Esto inyecta errores masivos a la matriz.

### La Solución (Clustering K-Means 1D)
El parámetro `--adaptive-rounding` (o `-ar`) enciende un algoritmo No Supervisado que:
1. Mira la nube de datos poblacional cruda.
2. Identifica dónde se desplazó el grupo (ej. se da cuenta de que todos se movieron a `0.38`).
3. Asigna inteligentemente a todo ese grupo su valor biológico real (`0.25`).

**💻 Ejemplo de Comando:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 --adaptive-rounding
```
* **Explicación bioinformática:** En lugar de redondear rígidamente, el programa analiza la distribución de la población completa en el SNP para descubrir los verdaderos grupos (clústeres) y corregir el sesgo de la máquina.

---

## 🚀 Nivel 5: El Modo "Dios" (Todo en Uno)

Al combinar todas estas funciones, obtienes el cálculo de dosis más robusto y estadísticamente avanzado posible en bioinformática de poliploides:

**💻 Ejemplo de Comando Maestro:**
```bash
java -jar biocenicana.jar allele-dosage \
    -v mis_datos_crudos.vcf \
    -p 10 \
    -md 15 \
    --adaptive-rounding \
    --impute bsdp-knn \
    -k 7
```

**Flujo interno del programa con este comando:**
1. **Filtra:** Elimina cualquier genotipo de dudosa calidad (`-md 15`).
2. **Corrige el Sesgo:** Usa Machine Learning para agrupar los datos sobrevivientes y asignar dosis biológicas exactas esquivando los errores del secuenciador (`--adaptive-rounding`).
3. **Rescata Datos:** A todos los genotipos vacíos o eliminados, los reconstruye copiando el ADN de sus 7 familiares más parecidos (`--impute bsdp-knn -k 7`).
