# Informe Técnico Metodológico: Suite BioCenicana
## Caracterización Avanzada de la Estructura Poblacional y Relación Genómica

### 1. Introducción
La suite **BioCenicana** es una plataforma bioinformática diseñada para el análisis robusto de la diversidad genética en poblaciones agrícolas complejas. Su arquitectura permite procesar datos de secuenciación masiva (NGS) representados en formato VCF, transformándolos en métricas interpretables mediante el uso de algoritmos de reducción de dimensionalidad, aprendizaje no supervisado (clustering) y estadística genómica poblacional. El enfoque principal radica en la gestión de organismos poliploides, donde la dosis alélica juega un papel crítico en la interpretación de la varianza genética.

---

### 2. Metodología de Procesamiento Bioinformático

#### 2.1. Representación de la Dosis Alélica
Dada una matriz de genotipos $G$ de dimensiones $n \times m$, donde $n$ es el número de muestras y $m$ el número de marcadores (SNPs), cada elemento $g_{ij}$ representa la dosis del alelo alternativo.
En organismos poliploides con nivel de ploidía $P$, la dosis se estima mediante:
$$d_{ij} = \text{round}\left( \frac{ad_{alt}}{ad_{ref} + ad_{alt}} \times P \right)$$
Donde $ad$ representa la profundidad de lectura (Allele Depth) para los alelos de referencia y alternativo.

#### 2.2. Normalización de Marcadores (Transformación de Patterson)
Para asegurar que los marcadores tengan el mismo peso estadístico independientemente de su frecuencia, se aplica una estandarización:
$$Z_{ij} = \frac{g_{ij} - E[g_j]}{\sqrt{Var(g_j)}} = \frac{g_{ij} - P \cdot p_j}{\sqrt{P \cdot p_j(1 - p_j)}}$$
Donde $p_j$ es la frecuencia del alelo alternativo en el marcador $j$. Esta normalización es fundamental para la correcta aplicación del análisis de componentes principales.

---

### 3. Reducción de Dimensionalidad (PCA)
El Análisis de Componentes Principales (PCA) se utiliza para proyectar la variación genética global en un espacio de baja dimensionalidad. Se implementa mediante la **Descomposición en Valores Singulares (SVD)** de la matriz $Z$:
$$Z = U \Sigma V^T$$
- **$U$ (Vectores Singulares Izquierdos):** Contiene los autovectores de la matriz de covarianza de las muestras.
- **$\Sigma$ (Valores Singulares):** Los cuadrados de estos valores ($\sigma_i^2$) corresponden a los autovalores ($\lambda_i$).
- **Eigenvalues ($\lambda_i$):** Representan la varianza capturada por cada componente.

La **Varianza Explicada Acumulada** se define como:
$$VE_{cum}(k) = \frac{\sum_{i=1}^k \lambda_i}{\sum_{j=1}^{\min(n,m)} \lambda_j}$$

---

### 4. Algoritmos de Agrupamiento y Estructura (Clustering)

#### 4.1. K-Means y Método del Codo
Se implementa K-Means minimizando la suma de cuadrados intra-cluster (WCSS):
$$WCSS = \sum_{k=1}^K \sum_{x \in C_k} ||x - \mu_k||^2$$
La determinación del número óptimo de poblaciones ($K$) se realiza mediante la búsqueda del "codo" en la función de $WCSS(K)$, identificando el punto donde la ganancia en explicación de varianza deja de ser significativa.

#### 4.2. Modelos de Mezcla Gaussiana (GMM) para Ancestría
Para modelar la ancestría compartida, se utiliza GMM con el algoritmo **Expectation-Maximization (EM)**. A diferencia de K-Means, GMM asigna una probabilidad de pertenencia (Ancestry Proportions) a cada muestra:
$$p(x_i) = \sum_{k=1}^K \pi_k \mathcal{N}(x_i | \mu_k, \Sigma_k)$$
Donde $\pi_k$ es el peso del componente y $\mathcal{N}$ es la distribución normal multivariada.

---

### 5. Análisis Discriminante (DAPC)
El DAPC (Discriminant Analysis of Principal Components) se utiliza para maximizar la separación entre los grupos predefinidos por el clustering inicial.
1. Se transforman los datos originales a componentes principales (PC space).
2. Se realiza un Análisis Discriminante Lineal (LDA) para encontrar las funciones discriminantes que maximizan el ratio de varianza entre-grupos ($B$) sobre la varianza intra-grupos ($W$):
$$\lambda = \max \frac{a^T B a}{a^T W a}$$
Los ejes LD resultantes son las direcciones de máxima diferenciación genética neta.

---

### 6. Relación Genómica: Matriz de Kinship (VanRaden)
Se implementa el método de **VanRaden (2008)** para estimar el parentesco genómico. La matriz $K$ se calcula como:
$$K = \frac{(G - P)(G - P)^T}{\sum 2 p_j (1 - p_j)}$$
Donde $P$ es una matriz que contiene las frecuencias alélicas medias por marcador. Esta matriz es esencial para corregir la estructura poblacional en estudios de asociación genómica (GWAS) y para la predicción genómica.

---

### 7. Agrupamiento Jerárquico (UPGMA)
La topología de las relaciones genéticas se visualiza mediante un dendrograma UPGMA (Unweighted Pair Group Method with Arithmetic Mean). La distancia entre dos grupos $A$ y $B$ se define como el promedio de todas las distancias entre pares de sus miembros:
$$d(A, B) = \frac{1}{|A| \cdot |B|} \sum_{a \in A} \sum_{b \in B} d(a, b)$$
Este método preserva las distancias ultramétricas y es ampliamente aceptado para la representación de distancias genéticas $F_{st}$ o Euclidianas.

---

**BioCenicana Technical Report**  
*Cenicaña Bioinformática*
