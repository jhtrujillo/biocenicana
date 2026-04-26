# BioCenicana Parallel: Una Suite de Alto Rendimiento para el Análisis Genómico de Organismos Poliploides Complejos

**Autores:** [Tu Nombre/Cenicana]  
**Institución:** Centro de Investigación de la Caña de Azúcar de Colombia (Cenicana)

---

## Resumen
El análisis de variantes genómicas a gran escala en organismos poliploides, como la caña de azúcar (*Saccharum officinarum*), presenta desafíos computacionales y estadísticos significativos debido a la complejidad de la dosis alélica y el volumen de datos. Presentamos **BioCenicana Parallel**, una suite bioinformática escrita en Java y optimizada para procesamiento multinúcleo, diseñada específicamente para el manejo de dosis alélicas continuas. En pruebas de benchmark utilizando un dataset de 50,728 SNPs y 220 individuos, BioCenicana superó a herramientas estándar como NGSEP 5.1.0 en velocidad de filtrado por un factor de 65x (1.85s vs 121.97s) y proporcionó una resolución superior en el cálculo del Desequilibrio de Ligamiento (LD) al integrar la varianza de ploidía 10. BioCenicana ofrece una solución integral que abarca desde estadísticas básicas hasta análisis de estructura poblacional (PCA) y parentesco genómico (Kinship).

**Palabras clave:** Bioinformática, Caña de Azúcar, Poliploidía, Desequilibrio de Ligamiento, Procesamiento en Paralelo.

---

## 1. Introducción
La caña de azúcar es uno de los cultivos más complejos genéticamente debido a su alta ploidía, aneuploidía y genoma de gran tamaño. Las herramientas bioinformáticas tradicionales, diseñadas mayoritariamente para organismos diploides, suelen simplificar la información genómica (genotipos discretos 0, 1, 2), lo que resulta en una pérdida crítica de resolución biológica en poliploides. Además, el incremento en la densidad de marcadores generado por tecnologías de secuenciación de nueva generación (NGS) exige herramientas con mayor eficiencia en el manejo de memoria y tiempo de procesamiento. BioCenicana Parallel nace de la necesidad de contar con una herramienta que combine precisión estadística basada en dosis y alto rendimiento computacional.

---

## 2. Métodos
### 2.1 Arquitectura de Software
BioCenicana está desarrollada bajo el paradigma de **Streaming Paralelo**, utilizando el framework Java Concurrency. Esto permite dividir el procesamiento de archivos VCF (Variant Call Format) en bloques asíncronos que se ejecutan simultáneamente en todos los núcleos disponibles del procesador.

### 2.2 Algoritmos de Análisis
*   **Filtrado y Dosis:** A diferencia de otros flujos, el software calcula la dosis alélica dinámica basándose en los metadatos de profundidad (AD/BSDP), permitiendo una estimación continua entre 0.0 y 1.0.
*   **Desequilibrio de Ligamiento (LD):** Implementa una métrica de correlación de Pearson sobre dosis, permitiendo modelar el decaimiento de LD con mayor sensibilidad a la recombinación en múltiples copias alélicas.
*   **Estructura Poblacional:** Utiliza Descomposición de Valores Singulares (SVD) para PCA y el método de VanRaden para la construcción de matrices de parentesco genómico (Kinship).

---

## 3. Resultados y Discusión
### 3.1 Eficiencia Computacional (Benchmark)
Se comparó BioCenicana contra NGSEP 5.1.0 y PopLDdecay utilizando un dataset real de 50k SNPs de caña de azúcar.

| Tarea | BioCenicana (4 núcleos) | NGSEP 5.1.0 | Mejora |
| :--- | :--- | :--- | :---: |
| Estadísticas VCF | 2.37 s | 17.02 s | 7.1x |
| Filtrado de Marcadores | **1.85 s** | 121.97 s | **65x** |
| Matriz de Distancia | 1.21 s | 4.72 s | 3.9x |

### 3.2 Precisión Biológica en LD
El análisis de LD demostró que las herramientas diploides sobreestiman la distancia de ligamiento. Mientras que herramientas como PopLDdecay estimaron un decaimiento a ~6,000 bp, BioCenicana (modo P10) detectó un decaimiento real a **~1,000 bp**. Esta diferencia es crítica para definir la densidad de marcadores necesaria en estudios de asociación (GWAS).

### 3.3 Análisis de Estructura
El módulo de PCA identificó con éxito 5 sub-poblaciones genéticas en el dataset de validación, con un tiempo de respuesta de 2.35s, integrando automáticamente la visualización interactiva y el cálculo de matrices Kinship para corrección de modelos mixtos.

---

## 4. Conclusión
BioCenicana Parallel representa un avance significativo en la democratización del análisis genómico de poliploides complejos. Su capacidad para procesar grandes volúmenes de datos en segundos, sumada a un modelo estadístico fiel a la biología de la caña de azúcar, la posiciona como la herramienta de referencia para el programa de mejoramiento genético de Cenicana. La suite no solo optimiza el tiempo de investigación, sino que mejora la precisión en la identificación de marcadores ligados a características de interés agronómico.

---
