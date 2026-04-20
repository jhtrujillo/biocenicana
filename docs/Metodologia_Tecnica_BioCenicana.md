# Manual de Referencia Científica: BioCenicana
## Metodología Avanzada de Análisis Poblacional y Estadística Genómica en Poliploides

### 1. Marco Teórico y Justificación
El análisis de la estructura poblacional es fundamental en genética de poblaciones para corregir el sesgo de estratificación en estudios de asociación genómica (GWAS) y para comprender la historia evolutiva de las variedades comerciales. En organismos poliploides, como la caña de azúcar, la complejidad aumenta debido a que cada locus puede presentar múltiples copias de un mismo alelo, lo que introduce una variación continua en la dosis alélica en lugar de las categorías discretas (0, 1, 2) propias de los diploides.

---

### 2. Algoritmos de Pre-procesamiento y Normalización

#### 2.1. Estimación de Dosis mediante Inferencia por Profundidad (AD)
La dosis alélica estimada ($\hat{d}_{i,j}$) para el individuo $i$ en el marcador $j$ se calcula bajo el supuesto de una distribución binomial de lecturas. Siendo $a_{i,j}$ el conteo del alelo alternativo y $n_{i,j}$ el conteo total ($ref + alt$), la dosis para una ploidía $P$ se define como:

$$\hat{d}_{i,j} = \arg\max_{k \in \{0..P\}} \binom{n_{i,j}}{a_{i,j}} \left(\frac{k}{P}\right)^{a_{i,j}} \left(1 - \frac{k}{P}\right)^{n_{i,j} - a_{i,j}}$$

En la práctica, BioCenicana aproxima esta verosimilitud mediante el redondeo del ratio de frecuencias observado, corregido por el nivel de ploidía del organismo.

#### 2.2. Estandarización de la Matriz Genómica $Z$
Sea $G$ la matriz de dosis cruda. Para el análisis de componentes principales, se requiere que la matriz de covarianza sea comparable. Aplicamos la normalización de **Patterson et al. (2006)**:

$$Z_{i,j} = \frac{g_{i,j} - \mu_j}{\sqrt{p_j (1 - p_j) (P)}}$$

Donde $\mu_j$ es la media observada y $p_j$ es la frecuencia del alelo minoritario (MAF). Esta transformación asegura que cada SNP tenga una varianza unitaria, evitando que los marcadores con mayor polimorfismo dominen artificialmente el análisis de varianza global.

---

### 3. Reducción de Dimensionalidad y Descomposición Espectral (PCA)

#### 3.1. Descomposición en Valores Singulares (SVD)
En lugar de calcular la matriz de covarianza $C = \frac{1}{m} ZZ^T$ (lo cual es computacionalmente costoso), BioCenicana utiliza SVD sobre la matriz $Z$:

$$Z = U \Sigma V^T$$

- **$U \in \mathbb{R}^{n \times n}$**: Los vectores singulares izquierdos, que representan las coordenadas de las muestras en el espacio latente.
- **$\Sigma \in \mathbb{R}^{n \times m}$**: Matriz diagonal de valores singulares $\sigma_i$.
- **Autovalores ($\lambda_i$):** Se derivan como $\lambda_i = \sigma_i^2 / (m-1)$.

#### 3.2. Significancia Estadística de los Componentes
La varianza explicada por cada componente se somete a la distribución de **Tracy-Widom** bajo la hipótesis nula de ausencia de estructura. BioCenicana permite visualizar este decaimiento mediante el *Scree Plot*, facilitando la selección de componentes que capturan variación biológica real frente al ruido estocástico.

---

### 4. Modelado Probabilístico de la Estructura (GMM)

El agrupamiento mediante **Gaussian Mixture Models (GMM)** permite una partición probabilística del espacio genético. El modelo se define como una combinación lineal de $K$ densidades gaussianas:

$$p(\mathbf{x}) = \sum_{k=1}^K \pi_k \mathcal{N}(\mathbf{x} | \boldsymbol{\mu}_k, \mathbf{\Sigma}_k)$$

Donde $\pi_k$ son los coeficientes de mezcla. La optimización se realiza mediante el algoritmo **Expectation-Maximization (EM)**:
1.  **E-Step:** Calcula la responsabilidad de cada componente $k$ para cada muestra $i$ ($q_{ik}$).
2.  **M-Step:** Actualiza los parámetros $\mu_k$ y $\Sigma_k$ maximizando la log-verosimilitud:
    $$\ln p(\mathbf{X} | \boldsymbol{\pi}, \boldsymbol{\mu}, \mathbf{\Sigma}) = \sum_{i=1}^n \ln \left( \sum_{k=1}^K \pi_k \mathcal{N}(\mathbf{x}_i | \boldsymbol{\mu}_k, \mathbf{\Sigma}_k) \right)$$

Este método permite estimar las **Ancestry Proportions** de forma análoga a los softwares *Admixture* o *Structure*, pero con una eficiencia computacional superior en el espacio de componentes.

---

### 5. Análisis Discriminante de Componentes Principales (DAPC)

Propuesto por **Jombart et al. (2010)**, el DAPC busca las combinaciones lineales de PCs que maximizan la separación entre grupos. El algoritmo resuelve el problema de autovalores para el ratio:

$$\max \frac{\mathbf{a}^T \mathbf{B} \mathbf{a}}{\mathbf{a}^T \mathbf{W} \mathbf{a}}$$

Donde:
- $\mathbf{B}$ es la matriz de dispersión **entre grupos** (Between-group scatter).
- $\mathbf{W}$ es la matriz de dispersión **dentro de los grupos** (Within-group scatter).

El resultado son los **Discriminantes Lineales (LDs)**, que representan los ejes donde las poblaciones están más diferenciadas genéticamente.

---

### 6. Matriz de Parentesco Genómico (VanRaden Kinship)

Para capturar la relación de identidad por descendencia (IBD), BioCenicana implementa el método de **VanRaden (2008)**, fundamental para el cálculo del BLUP genómico (gBLUP):

$$\mathbf{K} = \frac{\mathbf{M} \mathbf{M}^T}{\sum_{j=1}^m 2 p_j (1 - p_j)}$$

Aquí, $\mathbf{M}$ es la matriz de genotipos centrada ($G_{ij} - 2p_j$). El denominador es un factor de normalización basado en la heterocigosidad esperada bajo equilibrio de Hardy-Weinberg, lo que permite que los valores de la diagonal de $\mathbf{K}$ se aproximen a $1 + f$ (donde $f$ es el coeficiente de endogamia).

---

### 7. Análisis Filogenético y Distancia Genética

#### 7.1. Distancia de Cavalli-Sforza y Edwards
BioCenicana permite calcular distancias genéticas basadas en la geometría de la esfera, ideales para SNPs con deriva genética:
$$D_{chord} = \frac{2}{\pi} \sqrt{2(1 - \sum \sqrt{p_{1i} p_{2i}})}$$

#### 7.2. Algoritmo UPGMA
La reconstrucción jerárquica se realiza mediante un proceso iterativo de fusión de clusters que preserva la proporcionalidad del tiempo de divergencia, asumiendo una tasa de mutación constante (reloj molecular).

---

**Documento de Referencia Técnica**  
*División de Biometría y Bioinformática | Cenicaña*  
*Versión 1.1 - 2026*
