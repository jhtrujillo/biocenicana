# 🧬 Manual de Usuario: Allele Dosage Calculator (Biocenicana)

Este manual explica de forma clara y sencilla el contexto biológico, la matemática y el uso de las herramientas avanzadas de Inteligencia Artificial que incorpora `AlleleDosageCalculator`.

---

## 1. El Contexto Biológico: ¿Qué es la "Dosis Alélica"?

En organismos **diploides** (como los humanos), heredamos 2 copias de cada cromosoma (una de mamá, una de papá). Por lo tanto, si miramos una mutación (SNP), solo hay 3 posibilidades genéticas:
- `0.0`: Tienes 0 copias de la mutación.
- `0.5`: Tienes 1 copia de la mutación (Heterocigoto).
- `1.0`: Tienes 2 copias de la mutación.

Sin embargo, organismos como la **caña de azúcar, la papa o el trigo** son **poliploides**. La caña de azúcar moderna, por ejemplo, tiene ¡hasta 10 o 12 copias de cada cromosoma! 
Si la ploidía es 10, un individuo puede tener 0, 1, 2, 3... o 10 copias de esa mutación.

A la proporción de copias mutadas respecto al total se le llama **Dosis Alélica (Allele Dosage)**.
Matemáticamente, la calculamos usando las "lecturas" que hace la máquina de secuenciación:
`Dosis = Lecturas de Referencia / (Lecturas de Referencia + Lecturas Alternas)`

---

## 2. Fase 1: Filtro Estadístico de Profundidad (`--min-depth`)

### El Problema (Ruido Biológico)
Imagina que la máquina de secuenciación solo capturó **2 fragmentos de ADN** para un individuo en un SNP específico. Uno era Referencia, el otro Alterno. La matemática cruda diría: `1 / (1 + 1) = 0.5`. 
Pero si el organismo tiene 10 copias (ploidía 10), es imposible estar seguros con solo 2 lecturas de si su dosis real es `0.5`, `0.6` o `0.1`. Asignarle `0.5` basándose en solo 2 lecturas es **estadísticamente irresponsable**.

### La Solución
El parámetro `--min-depth` (o `-md`) descarta cualquier cálculo que no tenga un número mínimo de lecturas, marcándolo como "dato faltante" (`missing`) para que sea imputado de forma inteligente después.

**Ejemplo de uso:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 10 -md 10
```
*Traducción: "Calcula la dosis asumiendo ploidía 10, pero si un individuo tiene menos de 10 lecturas en total, descarta sus lecturas y márcalo como dato faltante".*

---

## 3. Fase 2: Imputación Inteligente KNN (`--impute knn`)

### El Problema (Imputación Clásica)
Si descartamos los datos dudosos (como vimos en la Fase 1), ¿qué ponemos en esos huecos vacíos?
El método clásico (`--impute mean` o `--impute mode`) rellena el hueco con el promedio poblacional.
- **Ejemplo:** Si tu población tiene cañas de la familia A (resistentes a un hongo) y cañas de la familia B (susceptibles), promediar la población entera para rellenar un hueco mezcla genes de familias diferentes y borra la estructura genética.

### La Solución (K-Nearest Neighbors - ML)
Activando `knn`, el programa usa un modelo de Inteligencia Artificial que funciona en dos pasos:
1. Analiza el genoma entero y construye un "árbol genealógico" virtual (Matriz de Distancias).
2. Cuando encuentra un hueco, busca a los **K individuos genéticamente más parecidos** (sus "vecinos/primos") y promedia solo las dosis de ellos.

**Ejemplo de uso:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 10 --impute knn -k 5
```
*Traducción: "Usa el método KNN para imputar huecos, basándote en los 5 parientes genéticos más cercanos de cada muestra".*

---

## 4. Fase 3: Redondeo Adaptativo (`--adaptive-rounding`)

### El Problema (Sesgo de Referencia)
En un mundo perfecto de un tetraploide (ploidía 4), los valores biológicos son fijos: `0.0, 0.25, 0.50, 0.75, 1.0`.
Pero los secuenciadores reales tienen un problema llamado **Reference Bias**: "leen" mejor el ADN original que el mutado.
- **Ejemplo:** A un individuo que biológicamente es `0.25`, la máquina lo reporta como `0.38` por el sesgo.
El método clásico matemático diría: *"0.38 está más cerca de 0.50 que de 0.25, le pondré 0.50"*. **¡Error genético grave!**

### La Solución (Clustering K-Means 1D)
El redondeo adaptativo ignora los cortes matemáticos fijos. En su lugar:
1. Mira la nube de datos de cómo salieron *todas* las muestras en ese SNP.
2. Descubre que "el grupo" se desplazó accidentalmente a `0.38`.
3. Se da cuenta de que `0.38` es el nuevo "centroide" para el nivel biológico `0.25`.
4. Asigna a todos los que cayeron en esa nube la dosis de `0.25`.

**Ejemplo de uso:**
```bash
java -jar biocenicana.jar allele-dosage -v archivo.vcf -p 4 --adaptive-rounding
```
*Traducción: "No uses cortes matemáticos rígidos. Permite que el algoritmo de Clustering agrupe a los individuos dinámicamente para corregir el sesgo de secuenciación".*

---

## 🚀 Poniendo todo junto (El Modo "Dios")

Puedes encadenar estas tres tecnologías para tener el cálculo de dosis más robusto, estadísticamente confiable y corregido que existe actualmente en la bioinformática de poliploides:

```bash
java -jar biocenicana.jar allele-dosage \
    -v mis_datos_crudos.vcf \
    -p 10 \
    -md 15 \
    --adaptive-rounding \
    --impute bsdp-knn \
    -k 7
```

**¿Qué está haciendo este comando paso a paso?**
1. Asume que trabajas con caña de azúcar (Ploidía 10).
2. **(Fase 1):** Elimina sin piedad a cualquier genotipo con menos de 15 lecturas (`-md 15`), porque no confiamos en ellos.
3. **(Fase 3):** Toma a los sobrevivientes y usa IA No Supervisada (`--adaptive-rounding`) para asignarles su dosis real, corrigiendo el sesgo de la máquina secuenciadora.
4. **(Fase 2):** A todos los genotipos que eliminamos en el paso 2, los reconstruye (`--impute bsdp-knn`) copiando el ADN de sus 7 familiares más parecidos (`-k 7`).
