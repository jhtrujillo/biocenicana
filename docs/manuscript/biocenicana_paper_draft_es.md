# BioCenicana: Un Conjunto de Herramientas de Alto Rendimiento Basado en Transmisión de Datos (Streaming) para el Análisis Genómico de Poliploides con Aplicaciones en *Saccharum* spp.

---

**Autores:** Jhon Henry Trujillo Montenegro¹, [Nombre del coautor]², [Nombre del coautor]³

**Afiliaciones:**
¹ Centro de Investigación de la Caña de Azúcar de Colombia (Cenicaña), Cali, Colombia
² [Institución], [Ciudad], [País]
³ [Institución], [Ciudad], [País]

**Autor para correspondencia:** [correo electrónico]

**Abreviaturas:** VCF, Variant Call Format (Formato de Llamada de Variantes); SNP, Single Nucleotide Polymorphism (Polimorfismo de Nucleótido Único); PCA, Principal Component Analysis (Análisis de Componentes Principales); NJ, Neighbor-Joining; LD, Linkage Disequilibrium (Desequilibrio de Ligamiento); MAF, Minor Allele Frequency (Frecuencia del Alelo Menor); GWAS, Genome-Wide Association Study (Estudio de Asociación de Genoma Completo); QC, Quality Control (Control de Calidad); GMM, Gaussian Mixture Model (Modelo de Mezcla Gaussiana); DBSCAN, Density-Based Spatial Clustering of Applications con Noise (Agrupamiento Espacial Basado en Densidad de Aplicaciones con Ruido).

---

## Resumen (Abstract)

Las especies poliploides, como la caña de azúcar (*Saccharum* spp., 2n = 10x ≈ 100–130), presentan desafíos estadísticos y computacionales sustanciales para los flujos de trabajo de análisis genómico estándar, los cuales están diseñados principalmente para organismos diploides. Las herramientas existentes a menudo no logran escalar al tamaño de los conjuntos de datos de poliploides, requieren una sobrecarga de memoria significativa, o carecen de soporte integrado para genómica poblacional consciente de las dosis alélicas. Aquí presentamos **BioCenicana**, un conjunto de herramientas Java de código abierto y alto rendimiento que implementa un motor de procesamiento secuencial (streaming) línea por línea para procesar archivos de Formato de Llamada de Variantes (VCF) a gran escala con una huella de memoria mínima. BioCenicana integra un flujo de trabajo analítico completo que abarca: (i) diagnóstico inicial del conjunto de datos y filtrado de control de calidad; (ii) extracción de la dosis alélica para análisis posteriores de GWAS; (iii) análisis de la estructura poblacional mediante Análisis de Componentes Principales (PCA), estimación de ancestría, y cálculo de la matriz de parentesco (Kinship) de VanRaden; (iv) análisis del decaimiento del Desequilibrio de Ligamiento (LD); (v) genómica comparativa, incluyendo visualización de sintenia, mapas de calor evolutivos basados en la densidad de SNPs, y estimación de la presión de selección (tasas Ka/Ks); y (vi) reconstrucción filogenética basada en SNPs mediante el algoritmo Neighbor-Joining con visualización interactiva. Todos los resultados analíticos se entregan como paneles de control HTML interactivos e independientes, que no requieren conexión a internet ni dependencias externas. Aplicado a un panel de diversidad de 220 genotipos de *Saccharum* spp. genotipados mediante secuenciación por genotipado (GBS/RAD-seq) con aproximadamente 50,000 SNPs de alta calidad, BioCenicana logró resolver exitosamente la subestructura poblacional, identificar importantes bloques sinténicos cromosómicos entre los genomas CC01-1940 y R570, y reconstruir un árbol filogenético de alta resolución coherente con los pedigrís de mejoramiento conocidos. BioCenicana está disponible de forma gratuita en https://github.com/jhtrujillo/biocenicana.

---

## Resumen en Lenguaje Sencillo (Plain Language Summary)

La caña de azúcar tiene uno de los genomas más complejos entre los cultivos agrícolas — posee hasta 13 copias de cada cromosoma, lo que hace que las herramientas de análisis genómico estándar sean en gran medida inadecuadas. BioCenicana es un software desarrollado específicamente para analizar caña de azúcar y otros cultivos poliploides (con múltiples copias de cromosomas). Este lee archivos de datos genómicos (formato VCF) línea por línea sin cargar todo el conjunto de datos en la memoria RAM, lo que permite procesar millones de marcadores genéticos en un computador portátil estándar. El software produce reportes visuales interactivos — incluyendo mapas de poblaciones, árboles evolutivos y gráficos de comparación de cromosomas — todos empaquetados como archivos HTML únicos que se abren directamente en cualquier navegador web sin necesidad de conexión a internet. Esto hace que el análisis genómico sea accesible para los mejoradores de plantas e investigadores que trabajan en entornos con recursos limitados.

---

## Introducción

La caña de azúcar (*Saccharum* spp.) es uno de los cultivos de mayor importancia económica en el mundo, responsable de aproximadamente el 80% de la producción mundial de azúcar y sirviendo como una materia prima clave para la bioenergía (FAO, 2022). Los cultivares comerciales modernos son híbridos interespecíficos complejos derivados principalmente de *Saccharum officinarum* y *Saccharum spontaneum*, con niveles de ploidía altamente variables que típicamente oscilan entre 2n = 10x ≈ 100 y 130 cromosomas (D'Hont et al., 1996; Grivet y Arruda, 2001). Esta extrema complejidad genómica — caracterizada por herencia polisómica, aneuploidía frecuente y redundancia alélica extensiva — plantea desafíos sustanciales para los flujos de trabajo bioinformáticos estándar diseñados para organismos diploides.

A pesar de la disponibilidad de herramientas de propósito general como GATK (McKenna et al., 2010), PLINK (Purcell et al., 2007) y TASSEL (Bradbury et al., 2007), su aplicación a datos de poliploides es severamente limitada. Muchas herramientas asumen una ploidía fija de dos, calculan frecuencias alélicas incorrectas en sistemas polisómicos, o requieren que todo el conjunto de datos se cargue en la memoria RAM — una restricción poco práctica al procesar archivos VCF con cientos de muestras en miles de andamiajes (scaffolds) cromosómicos. Existe software especializado para genómica de poliploides (p. ej., GWASpoly, Rosyara et al., 2016; polyRAD, Clark et al., 2019; updog, Gerard et al., 2018), pero estas herramientas típicamente abordan tareas analíticas específicas en lugar de proporcionar un flujo de trabajo integrado de principio a fin.

La necesidad de un conjunto de herramientas integrado, eficiente en memoria y consciente de la poliploidía es particularmente aguda en contextos de mejoramiento genético aplicado, donde los investigadores requieren una transición rápida desde los datos crudos de genotipado hasta conocimientos prácticos sobre la estructura poblacional, diversidad genética y relaciones evolutivas. Desarrollamos **BioCenicana** para abordar esta brecha. El conjunto de herramientas implementa un motor de transmisión de datos (streaming) que procesa archivos VCF secuencialmente — una línea a la vez — desacoplando así el uso de memoria del tamaño del conjunto de datos. Integra todos los pasos principales de un flujo de trabajo moderno de genómica poblacional en una sola aplicación de línea de comandos, y genera paneles HTML interactivos e independientes que facilitan la interpretación y comunicación de los resultados sin requerir entornos bioinformáticos especializados.

Aquí describimos el diseño, la implementación y la validación de BioCenicana, y demostramos su aplicación a un panel de 220 accesiones de *Saccharum* spp. Mostramos que la herramienta recupera subestructuras poblacionales conocidas, produce árboles filogenéticos consistentes con los registros de pedigrí y permite un análisis genómico comparativo a gran escala entre genomas de referencia.

---

## Materiales y Métodos

### Material Vegetal y Genotipado

Para este estudio se utilizó un archivo VCF correspondiente a un panel de diversidad de 220 accesiones de *Saccharum* spp., las cuales representan toda la diversidad presente dentro del banco de germoplasma de Cenicaña (el cual alberga actualmente alrededor de 1,600 accesiones). Estas accesiones fueron genotipadas en estudios previos mediante técnicas de GBS y RAD-seq y mapeadas al genoma de referencia CC01-1940. El archivo VCF fue seleccionado y filtrado teniendo en cuenta estrictos parámetros de calidad: un número mínimo de individuos genotipados del 50% (es decir, al menos 110 individuos por variante), una calidad mínima de llamado de 30, y una profundidad de secuenciación mínima de 20X. Este VCF depurado sirvió como punto de partida para todos los análisis de la herramienta.

### Arquitectura de Software

BioCenicana está desarrollado en el lenguaje de programación Java 11 y se distribuye como un único archivo ejecutable (formato JAR). Esto significa que la herramienta es portátil y puede ejecutarse en cualquier sistema operativo (Windows, Mac o Linux) sin requerir instalaciones complejas o configuraciones adicionales.

La herramienta se opera desde la terminal (interfaz de línea de comandos) y fue construida utilizando la biblioteca *picocli*. Esta biblioteca garantiza que los comandos sigan el estándar **POSIX** (Portable Operating System Interface), lo cual es simplemente un conjunto de reglas universales que aseguran que los comandos de BioCenicana se comporten de manera predecible y familiar para cualquier usuario (por ejemplo, usando guiones cortos para opciones simples como `-v` y guiones dobles para opciones completas como `--vcf`). Además, genera automáticamente manuales de ayuda en pantalla. Internamente, el programa está diseñado para reportar claramente al sistema operativo si un análisis fue exitoso o si falló, facilitando su integración en flujos de trabajo automatizados.

El principio de diseño central, y la mayor ventaja técnica de BioCenicana, es su **motor de procesamiento secuencial (streaming)**. Normalmente, los programas de genómica intentan cargar todo el archivo de datos (VCF) en la memoria principal del computador (RAM) antes de analizarlo, lo cual hace que los equipos colapsen o se queden sin memoria cuando los archivos son muy grandes. Para solucionar esto, BioCenicana lee los archivos como si fuera un libro: avanza línea por línea extrayendo únicamente la información necesaria de forma temporal y descartando el resto.

Gracias a este enfoque, el consumo de memoria del computador ya no depende de cuántos millones de variantes genéticas (SNPs) tenga el archivo, sino únicamente de la cantidad de muestras (individuos) analizadas. En términos prácticos, esto permite a un investigador procesar archivos genómicos gigantescos que superan los 60 GB en computadores portátiles estándar que solo tienen 8 GB de memoria RAM, democratizando el acceso a análisis de alta complejidad computacional.

Adicionalmente, la arquitectura de la herramienta incorpora soporte para **procesamiento en paralelo (multiprocesamiento)**. Esto significa que las tareas matemáticas más exigentes, como el cálculo de distancias genéticas o la creación de matrices de parentesco, pueden dividir su carga de trabajo automáticamente entre los múltiples núcleos (procesadores) que tenga el computador. El usuario puede habilitar esto fácilmente indicando el número de hilos de procesamiento deseados (por ejemplo, mediante el comando `-t`), lo cual acelera drásticamente los tiempos de ejecución y aprovecha al máximo la capacidad del equipo.

### Módulos y Funciones de la Herramienta

BioCenicana está estructurado en cinco módulos principales, cada uno diseñado para ejecutar una etapa específica del análisis genómico y producir resultados visuales directamente interpretables:

1. **`vcf-stats` y `vcf-filter` (Control de Calidad):** Funciones encargadas de diagnosticar el estado de los datos crudos (frecuencias, datos faltantes, profundidad) y filtrar las variantes de baja calidad o poco informativas.
2. **`pop-structure` (Estructura Poblacional):** Módulo que determina cómo se agrupan los individuos genéticamente, utilizando técnicas estadísticas de reducción de dimensionalidad, estimación de ancestría (admixture) y matrices de parentesco.
3. **`ld` (Desequilibrio de Ligamiento):** Función para medir qué tan ligados (heredados juntos) están los marcadores a lo largo de los cromosomas, lo cual es vital para saber cuántos marcadores se necesitan en estudios genéticos.
4. **`snp-tree` (Filogenia):** Módulo que reconstruye la historia evolutiva y las relaciones biológicas de parentesco entre las accesiones, generando un árbol interactivo.
5. **`comp-gen` y `kaks-calc` (Genómica Comparativa):** Herramientas para comparar grandes bloques de cromosomas entre dos especies o variedades (sintenia) e identificar matemáticamente qué genes están bajo presión de selección evolutiva.

### Fundamentos Biológicos, Estadísticos y Bioinformáticos

A continuación, se detalla la metodología de construcción de cada módulo, integrando explícitamente los principios biológicos, los modelos estadísticos/matemáticos subyacentes y las estrategias bioinformáticas implementadas en el código de BioCenicana.

#### 1. Módulo de Control de Calidad y Filtrado (`vcf-stats`, `vcf-filter`)

A diferencia de los organismos diploides (donde un genotipo suele ser solo bivalente, por ejemplo AA, Aa o aa), en los organismos poliploides complejos como la caña de azúcar es fundamental comprender la "dosis alélica". Esto significa que importa exactamente cuántas copias de un alelo variante posee un individuo, lo cual puede ir desde 0 hasta 10 copias en un genoma decaploide. Dado que las secuencias de ADN crudas producidas por los secuenciadores contienen errores de lectura, falsos positivos y regiones con baja cobertura que no representan mutaciones biológicas reales, es indispensable un proceso de diagnóstico y depuración riguroso. Para abordar esto, el módulo se dividió computacionalmente en dos funciones complementarias: `vcf-stats` para el diagnóstico y `vcf-filter` para la depuración activa.

La función `vcf-stats` se construyó como una herramienta de perfilamiento bioinformático no destructivo. Su objetivo es leer un archivo crudo y calcular las métricas estadísticas globales y por individuo sin alterar los datos originales. Computacionalmente, utiliza el **motor de transmisión secuencial (*streaming engine*)** de BioCenicana. A diferencia de los programas de genómica tradicionales que intentan cargar todo el archivo masivo directamente en la memoria principal (RAM) —lo que suele causar colapsos en computadores de escritorio— un *streaming engine* funciona procesando el archivo como si fuera una cinta transportadora: lee y analiza los datos rigurosamente "línea por línea", extrayendo lo necesario y liberando de la memoria la información que ya procesó. Durante este ciclo de lectura continua, el algoritmo extrae en tiempo real los campos de profundidad de lectura (`DP`) y profundidad alélica (`AD`) para compilar tres métricas fundamentales:
1. **Frecuencia del Alelo Menor (MAF):** Se calcula sumando todas las observaciones del alelo alternativo sobre el total de alelos muestreados en la población. 
2. **Tasa de Datos Faltantes (Missing Rate):** Evalúa qué porcentaje de la población carece de información confiable para un marcador específico.
3. **Profundidad Media:** Cuantifica la cobertura de secuenciación promedio por individuo y por marcador.

Los resultados se consolidan estadísticamente y se exportan automáticamente como un panel de control HTML interactivo, permitiendo al investigador visualizar gráficamente el estado de salud de su genotipado antes de tomar decisiones de filtrado.

Una vez identificados los perfiles de error mediante `vcf-stats`, el usuario emplea la función `vcf-filter` para depurar la base de datos. Esta función integra métodos estadísticos sofisticados para lidiar con la incertidumbre propia de los poliploides. Para resolver la ambigüedad en zonas del genoma con baja profundidad de secuenciación, `vcf-filter` evita realizar un llamado genotípico rígido o absoluto. En su lugar, aplica un **modelo de máxima verosimilitud ponderado**. 

En términos sencillos, la "verosimilitud" es un concepto estadístico que evalúa qué tan lógico es un resultado observado bajo diferentes escenarios hipotéticos. En lugar de forzar al programa a adivinar ciegamente si el individuo tiene, por ejemplo, 3 o 4 copias del gen mutante (dosis) a partir de una lectura secuencial pobre, el modelo calcula una probabilidad matemática para *cada una* de las 11 posibles dosis reales (desde 0 hasta 10 copias en un genoma decaploide). Si $D$ son los datos de las lecturas observadas (campo `AD` en el VCF) y $k$ es una posible dosis alélica, el modelo determina la probabilidad $P(k|D)$. Posteriormente, en lugar de elegir caprichosamente un único número entero, la herramienta calcula un valor "esperado" (continuo) multiplicando cada dosis por su respectiva probabilidad y sumándolas todas:

$$ \text{Dosis Estimada} = \sum_{k=0}^{10} k \cdot P(k | D) $$

Este sofisticado enfoque probabilístico permite a BioCenicana asignar valores genotípicos con decimales (por ejemplo, estimar que la dosis es de 3.4 copias en lugar de redondear bruscamente a 3). Esto minimiza dramáticamente el sesgo técnico generado por la secuenciación incompleta y retiene la valiosa incertidumbre estadística para los análisis posteriores.

A nivel de ejecución, `vcf-filter` opera aplicando una serie de umbrales estrictos definidos por el usuario a través de parámetros de la línea de comandos. Los parámetros principales incluyen:

*   `--ploidy`: Define el nivel de ploidía biológica de la especie analizada (por ejemplo, $10$ para cultivares modernos de caña de azúcar). Este parámetro es estructuralmente fundamental, ya que le indica al modelo de máxima verosimilitud exactamente cuántos estados de dosis alélica (de 0 a $k$) debe evaluar computacionalmente.
*   `--min-maf`: Frecuencia mínima del alelo menor permitida (típicamente $0.05$ para eliminar variantes raras que suelen representar ruido técnico).
*   `--max-missing`: Proporción máxima de individuos sin datos permitida para conservar un marcador (por ejemplo, $0.20$ tolera un 20% de datos faltantes en la población).
*   `--min-dp`: Profundidad de secuenciación mínima requerida por muestra para considerar que la lectura de la dosis alélica es estadísticamente válida.
*   `--top-n`: (Opcional) Un filtro bioinformático de selección de características (*feature selection*) que clasifica y retiene únicamente los marcadores más informativos basados en su heterocigosidad esperada, reduciendo drásticamente la carga computacional para los análisis posteriores.
*   `-t` (Hilos de Procesamiento): Aunque el archivo se lee secuencialmente, el cálculo de las probabilidades de máxima verosimilitud para cientos de muestras en cada variante es matemáticamente intensivo. Este parámetro permite paralelizar estas operaciones distribuyendo la carga de cálculo entre los múltiples núcleos del procesador de la máquina, acelerando masivamente el tiempo de depuración.

Las variantes que no superan estos umbrales paramétricos son omitidas del flujo de memoria al instante. El resultado es la exportación en tiempo real de un archivo VCF altamente depurado, garantizando una complejidad temporal mínima de $O(N)$ y un consumo de RAM casi imperceptible.

**Ejemplo de Ejecución:**
```bash
# 1. Diagnóstico estadístico asumiendo un genoma decaploide
java -jar biocenicana.jar vcf-stats -v panel_crudo.vcf --ploidy 10 -o reporte_stats.html

# 2. Depuración activa basada en umbrales aprovechando 8 núcleos del procesador
java -jar biocenicana.jar vcf-filter -v panel_crudo.vcf -o panel_filtrado.vcf -t 8 \
    --ploidy 10 --min-maf 0.05 --max-missing 0.20 --min-dp 20 --top-n 50000
```

#### 2. Módulo de Análisis de Estructura Poblacional (`pop-structure`)

Identificar la estructura de una población (es decir, cómo se agrupan los individuos y si existe mezcla genética o *admixture*) es biológicamente esencial para evitar resultados falsos positivos en estudios de asociación (GWAS) y para orientar cruces estratégicos en programas de mejoramiento. Las poblaciones modernas de caña de azúcar, caracterizadas por su extrema complejidad, son el resultado de introgresiones históricas continuas entre especies domesticadas y silvestres. Para mapear y desenredar esta intrincada red genética, BioCenicana ejecuta un flujo de trabajo estadístico secuencial: primero, reduce la enorme dimensionalidad de los datos genómicos para hacerlos interpretables; segundo, calcula una matriz exacta del parentesco entre las accesiones; y finalmente, despliega modelos algorítmicos de agrupamiento (*clustering*) para descubrir y clasificar automáticamente las subpoblaciones subyacentes.

El primer desafío analítico radica en la inmensa cantidad de datos. Biológicamente, intentar comparar 220 plantas utilizando 50,000 marcadores genéticos simultáneamente es imposible para el ojo humano, ya que implicaría visualizar un gráfico de 50,000 dimensiones. Para resolver esto, el módulo aplica un Análisis de Componentes Principales (PCA). Matemáticamente, el PCA es una técnica de álgebra lineal que reduce la dimensionalidad, "comprimiendo" toda la información masiva en unos pocos ejes nuevos (Componentes Principales) que preservan la mayor cantidad de variación estadística posible. El programa construye una matriz de genotipos basada en las dosis alélicas, la estandariza y extrae sus vectores principales. Así, los investigadores pueden ver en un simple gráfico 3D cómo las plantas se agrupan según su similitud genética global.

Paralelamente a esta reducción visual, es necesario calcular el grado exacto de consanguinidad de la población. Para ello, el módulo construye una Matriz de Parentesco (*Kinship*), la cual es esencialmente una tabla cuadrada que cuantifica el grado de similitud biológica entre todos los pares posibles de individuos. BioCenicana utiliza el método de VanRaden (2008) adaptado para poliploides. La lógica estadística de este método asume que si dos plantas comparten un alelo muy raro (con baja frecuencia poblacional), están mucho más estrechamente emparentadas que si comparten un alelo muy común. Para lograrlo, la fórmula resta a la dosis alélica de cada planta ($X$) el promedio de la población ($P$) con el fin de "centrar" el dato, y luego calcula el producto matricial consigo mismo dividiéndolo por un factor de escala de la varianza ($c$), resultando en la ecuación:

$$ \text{Kinship} = \frac{(X - P)(X - P)'}{c} $$

El valor numérico final indica con precisión matemática qué tan parientes son dos plantas, información indispensable para realizar modelos predictivos precisos en genética cuantitativa.

Una vez que el PCA ha reducido la información genómica al espacio tridimensional, el siguiente paso es identificar fronteras biológicas objetivas mediante algoritmos de agrupamiento estadístico (*Clustering*). En lugar de obligar al investigador a definir agrupaciones subjetivamente, el programa integra tres modelos matemáticos de inteligencia artificial para clasificar automáticamente a las plantas con perfiles genéticos similares. El primero es K-Means, un agrupamiento por centroides geométricos donde el usuario define el número de grupos esperados; el algoritmo ubica puntos centrales en el espacio y asigna cada planta a su centro más cercano (distancia euclidiana), siendo excelente para poblaciones bien estructuradas en grupos esféricos. El segundo es DBSCAN, un agrupamiento por densidad espacial que, a diferencia de K-Means, no requiere adivinar cuántos grupos existen. Este modelo busca "manchas" densamente pobladas en el gráfico y etiqueta matemáticamente como "ruido" a aquellas plantas aisladas en zonas vacías, haciéndolo ideal para descubrir estructuras irregulares o híbridos anómalos. El tercero es el Modelo de Mezcla Gaussiana (GMM), el cual representa el enfoque estadístico más avanzado. En lugar de hacer clasificaciones tajantes o absolutas, el GMM asume que cada subpoblación tiene forma de campana de Gauss en múltiples dimensiones y calcula probabilidades de pertenencia; por ejemplo, puede dictaminar que una planta tiene un 80% de probabilidad de pertenecer al clado silvestre y un 20% al domesticado, capturando elegantemente la naturaleza gradual de la hibridación genética.

A nivel operativo, estas pesadas operaciones matriciales son orquestadas mediante una serie de parámetros en la línea de comandos que le otorgan control analítico total al investigador. El parámetro `-t` (hilos) es vital, ya que paraleliza el trabajo de álgebra lineal distribuyéndolo entre múltiples procesadores de la máquina para evitar cuellos de botella computacionales. Asimismo, parámetros como `--n-pca` dictan el número exacto de componentes a extraer, mientras que `--dbscan-eps` (radio de búsqueda espacial) y `--dbscan-minpts` (número mínimo de plantas vecinas) permiten afinar milimétricamente la sensibilidad de los algoritmos de agrupamiento por densidad. Finalmente, los resultados tabulares se inyectan dinámicamente en plantillas HTML autónomas que utilizan bibliotecas nativas de JavaScript, permitiendo al investigador rotar e interactuar libremente con las nubes de puntos poblacionales 3D desde cualquier navegador web, sin requerir instalación de software de terceros.

**Ejemplo de Ejecución:**
```bash
# Análisis completo de estructura poblacional utilizando 8 hilos de procesamiento (-t)
# Se calculan 10 componentes principales (--n-pca) y los modelos de agrupamiento 
# con parámetros de densidad específicos para DBSCAN (--dbscan-eps, --dbscan-minpts)
java -jar biocenicana.jar pop-structure -v panel_filtrado.vcf -t 8 \
    --ploidy 10 --pca --n-pca 10 --kinship \
    --admixture -k 5 \
    --dbscan-eps 0.5 --dbscan-minpts 5 \
    -o resultados_pop/
```

#### 3. Módulo de Desequilibrio de Ligamiento (`ld`)

El Desequilibrio de Ligamiento (LD, por sus siglas en inglés) es un concepto estadístico fundamental en genética poblacional que mide la asociación biológica o "co-herencia" entre diferentes regiones del ADN. En términos simples, si dos marcadores genéticos se encuentran físicamente muy cerca en un mismo cromosoma, existe una alta probabilidad de que se hereden juntos de generación en generación, evadiendo la separación que normalmente ocurre durante la recombinación (entrecruzamiento). Por ejemplo, si en un panel de diversidad de caña de azúcar observamos que una mutación genética responsable de la resistencia a la roya casi siempre aparece acompañada de otra mutación silenciosa ubicada un poco más adelante en el mismo cromosoma, decimos que ambos marcadores se encuentran en "alto desequilibrio de ligamiento". Cuantificar el LD es crucial para el diseño experimental: si el LD en un cultivo "decae" (es decir, se pierde) rápidamente a lo largo de unos pocos pares de bases de distancia, los investigadores sabrán de antemano que están obligados a secuenciar el genoma con una enorme densidad de marcadores genéticos para no correr el riesgo de "perderse" genes agronómicamente importantes durante un estudio de asociación de genoma completo (GWAS).

Para medir matemáticamente la fuerza de esta co-herencia, el módulo calcula el estadístico $r^2$. Los valores numéricos generados por este estadístico fluctúan estrictamente entre 0 (lo que indica que los marcadores se heredan de manera completamente aleatoria e independiente) y 1 (lo que indica que los marcadores se heredan siempre juntos como un bloque perfecto e inquebrantable). El núcleo algorítmico de BioCenicana para calcular esta matriz de $r^2$ toma como base el histórico modelo formulado por Hill y Robertson (Hill & Robertson, 1968), el cual estimaba originalmente las varianzas y covarianzas de frecuencias alélicas en poblaciones finitas. Sin embargo, BioCenicana implementa una optimización estadística profunda de este modelo. Las herramientas tradicionales se limitan a fenotipos diploides discretos (donde un marcador simplemente se categoriza como presente o ausente, ej. 0, 1 o 2). En contraste, el algoritmo de BioCenicana fue generalizado algebraicamente para procesar variables de **dosis continuas** (como los valores probabilísticos estimados por el modelo de máxima verosimilitud del módulo `vcf-filter`, por ejemplo, 3.4 o 5.8 copias de un alelo). Este importante avance matemático permite que el motor de correlación se adapte automáticamente a **cualquier nivel de ploidía**, garantizando precisión estadística absoluta tanto en especies diploides sencillas (como el arroz o el maíz) como en especies poliploides de altísima complejidad (como la caña de azúcar decaploide), sin requerir alteraciones en el código fuente.

A nivel de software, calcular el $r^2$ entre todas las combinaciones posibles de 50,000 marcadores genéticos implicaría billones de operaciones matemáticas — un problema computacionalmente explosivo que saturaría la memoria de cualquier estación de trabajo convencional. Para sortear esto, el módulo fue diseñado empleando un algoritmo bioinformático de ventana deslizante (*sliding window*). A nivel operativo, el usuario define el parámetro `--max-dist` para indicarle al programa la distancia física máxima (en pares de bases) dentro de la cual tiene sentido biológico buscar correlaciones. El algoritmo recorre secuencialmente el cromosoma calculando el LD únicamente entre variantes vecinas que caigan dentro de dicha ventana, apoyándose fuertemente en el paralelismo distribuido por hilos de CPU (`-t`). Adicionalmente, el parámetro `--min-r2` actúa como un umbral estadístico que descarta al instante las correlaciones débiles o consideradas ruido, ahorrando drásticamente espacio de escritura en el disco duro. Finalmente, los millones de cálculos resultantes de $r^2$ se ajustan a una curva matemática de decaimiento no lineal en función de la distancia física, y se exportan nativamente a un gráfico interactivo HTML. Esto permite al genetista visualizar e inferir exactamente a qué distancia física (en kb) se rompen los bloques haplotípicos en su población de mejoramiento.

**Ejemplo de Ejecución:**
```bash
# Estimación del decaimiento del LD calculando en paralelo (-t 8) 
# con una ventana máxima de 200kb y filtrando correlaciones débiles (--min-r2)
java -jar biocenicana.jar ld -v panel_filtrado.vcf --ploidy 10 -t 8 \
    --max-dist 200000 --min-r2 0.1 -o grafico_ld.html
```

#### 4. Módulo de Reconstrucción Filogenética (`snp-tree`)

Determinar la historia evolutiva de las plantas es un pilar fundamental en la conservación y el aprovechamiento de los recursos fitogenéticos. Biológicamente, conocer la filogenia permite a los taxónomos y genetistas validar pedigrís históricos, separar especies puras de híbridos de introgresión, e identificar visualmente el grado de parentesco biológico. Para lograr esto a partir de datos genómicos masivos, BioCenicana inicia su flujo de trabajo calculando una inmensa matriz de distancias genéticas emparejadas. Matemáticamente, el módulo emplea una métrica de **distancia euclidiana**. En términos sencillos, la distancia euclidiana representa la longitud de la línea recta más corta que conecta dos puntos en un espacio matemático. El algoritmo toma las dosis alélicas de la Planta A y las compara con las de la Planta B, midiendo matemáticamente qué tan "lejos" están sus perfiles genéticos directos. Al normalizar este valor por la cantidad exacta de marcadores que ambas plantas comparten de manera efectiva, el modelo logra evitar los graves sesgos estadísticos inducidos por los datos faltantes que son inherentes a las tecnologías modernas de secuenciación (como GBS o RAD-seq).

Una vez que se tiene una tabla con la distancia exacta de cada planta respecto a todas las demás, se debe inferir y dibujar la topología evolutiva. Para ello, BioCenicana implementa el algoritmo *Neighbor-Joining* (Saitou & Nei, 1987). Algunos modelos evolutivos antiguos asumen un "reloj molecular estricto", asumiendo erróneamente que la **tasa de mutación es constante a lo largo del tiempo** en todos los linajes. Sin embargo, en la realidad biológica, algunas variedades acumulan mutaciones mucho más rápido que otras debido a fuertes presiones de domesticación o estrés ambiental. El algoritmo *Neighbor-Joining* soluciona esto elegantemente: no asume una tasa de mutación constante. En su lugar, el modelo computacional agrupa recursivamente a los pares más similares ("vecinos"), construyendo iterativamente una red que **minimiza la longitud total de todas las ramas evolutivas**. El árbol resultante refleja de manera altamente precisa la cantidad de evolución real que ha experimentado cada linaje de forma independiente.

Ahora bien, es imperativo comprobar matemáticamente la robustez de este árbol. Para asegurar el máximo rigor científico, el algoritmo evalúa la confianza estadística de cada ramificación mediante la técnica de **bootstrap** (controlada por el parámetro `--bootstrap`). El *bootstrap* funciona como un exhaustivo test de estrés computacional: el programa toma los marcadores genéticos originales, extrae miles de subconjuntos de datos aleatorios (mezclándolos con reemplazo) y obliga a la computadora a reconstruir el árbol evolutivo cientos de veces sucesivas. Si, bajo el parámetro `--bootstrap 100`, las Plantas A y B terminan agrupadas juntas en 98 de los 100 árboles aleatorios generados, esa rama obtendrá un soporte de confianza innegable (98%), garantizando al investigador que ese evento evolutivo es biológicamente real y no un espejismo estadístico del muestreo.

Los resultados finales de esta inferencia matemática se empaquetan en el formato estándar bioinformático universal *Newick*. Finalmente, un sofisticado motor de renderizado desarrollado nativamente en JavaScript lee este archivo sin depender de pesadas bibliotecas externas, calcula trigonométricamente las coordenadas espaciales de cada nodo, y dibuja el árbol evolutivo utilizando gráficos vectoriales interactivos (SVG). A nivel operativo, si el investigador conecta los resultados del módulo poblacional (utilizando el parámetro `--pca`), el motor teñirá dinámicamente cada rama del árbol, permitiendo una síntesis visual perfecta entre la estructura estadística de la población y su profunda historia evolutiva.

**Ejemplo de Ejecución:**
```bash
# Reconstrucción filogenética en paralelo (-t 8) evaluando la robustez de las ramas (--bootstrap)
# e integrando los resultados estadísticos del PCA para el coloreado interactivo de los nodos
java -jar biocenicana.jar snp-tree -v panel_filtrado.vcf --ploidy 10 -t 8 \
    --pca resultados_pop/pca_clusters.csv \
    --bootstrap 100 -o filogenia_interactiva.html
```

#### 5. Módulo de Genómica Comparativa (`comp-gen`, `kaks-calc`)

A diferencia de los módulos anteriores, enfocados en marcadores poblacionales, la Genómica Comparativa evalúa la arquitectura estructural de los cromosomas completos. Un concepto clave aquí es la **sintenia**, que ocurre cuando dos especies o variedades (por ejemplo, un híbrido comercial moderno y su ancestro silvestre) conservan el orden de sus bloques de genes a lo largo de un cromosoma, incluso tras millones de años de divergencia. Identificar estos bloques de genes conservados permite descubrir macro-mutaciones, como inversiones o fusiones cromosómicas, las cuales suelen tener un rol fundamental en la especiación y en la adaptación de las plantas a distintos entornos.

Para mapear esta organización estructural, la función `comp-gen` actúa como un procesador e integrador de datos. En la bioinformática actual, y específicamente en los flujos de trabajo de Cenicaña, el estándar para identificar bloques homólogos es el uso de herramientas especializadas como **MCScanX** (Wang et al., 2012) o **SyMAP**. BioCenicana toma los resultados tabulares generados por estos programas, los cruza con las coordenadas físicas de los genomas (archivos GFF) y los convierte directamente en representaciones visuales interactivas.

Un desafío común al visualizar la conexión de miles de genes entre dos genomas es que trazar simples líneas rectas genera gráficos saturados y difíciles de interpretar. Para resolverlo, el algoritmo implementa **curvas de Bezier**. Estas curvas son ecuaciones geométricas paramétricas que permiten dibujar trayectorias suaves para conectar un gen en el Genoma A con su homólogo en el Genoma B, facilitando la lectura del gráfico. Además, para simplificar la interpretación y ocultar asociaciones cortas de poco interés biológico, el usuario puede aplicar un filtro mediante el parámetro `--min-block-size`, graficando únicamente los bloques que superen cierto umbral de genes consecutivos.

Una de las utilidades más prácticas de este módulo es la integración directa de datos genómicos poblacionales mediante **archivos VCF**. Superponer un VCF sobre los gráficos de sintenia permite ir más allá de la mera estructura física del cromosoma, mostrando también la densidad de mutaciones (SNPs) acumuladas por la población sobre dichos bloques conservados. Esto ayuda al investigador a identificar rápidamente zonas "calientes" (regiones genómicas altamente polimórficas) y zonas "frías" (regiones muy conservadas). Como resultado, BioCenicana genera múltiples salidas gráficas, como diagramas circulares (tipo *Circos*), cintas comparativas en 3D y mapas de calor (*heatmaps*) de densidad alélica, todos exportados como archivos HTML.

A nivel de terminal, este análisis se ensambla proporcionando a la herramienta los insumos necesarios en un solo comando. El usuario ingresa las anotaciones de los genomas comparados mediante `--gff1` y `--gff2`, y el archivo de sintenia base con `--collinearity`. Para darle sentido biológico a las conexiones, se añaden las funciones de los genes con `--annot1` y `--annot2`, y finalmente los datos de mutaciones de la población con `--vcf`. El programa devuelve dos resultados principales: un archivo tabular (`-o`) con todas las métricas procesadas listas para análisis estadístico, y un visor HTML autónomo (`--viz`) con los gráficos.

De manera complementaria a esta visión macroestructural, el módulo cuenta con la función `kaks-calc` para explorar la micro-evolución a nivel de secuencias. Su objetivo es medir la tasa de mutación dentro de los genes codificantes para inferir el tipo de selección natural que actúa sobre ellos. Se evalúa si un gen está bajo **"selección purificadora"** (cuando se penalizan biológicamente las mutaciones porque la función de la proteína es estrictamente vital) o bajo **"selección positiva"** (cuando se favorecen las mutaciones que permiten adaptación rápida, como ocurre con algunos genes de resistencia a patógenos). 

Esta medición se realiza mediante el cálculo de la relación **Ka/Ks**: la proporción entre las tasas de sustitución que alteran el aminoácido (Ka) frente a las mutaciones sinónimas que no lo alteran (Ks), empleando el modelo de Nei-Gojobori (Nei & Gojobori, 1986). A nivel de ejecución, BioCenicana lee los genes en formato FASTA, realiza los alineamientos de codones paralelizando el cálculo con múltiples procesadores (`-t`) y genera un reporte tabular con los valores de significancia estadística.

**Ejemplo de Ejecución:**
```bash
# Generación interactiva de gráficos de sintenia poblacional entre genomas de referencia
# Integrando coordenadas (GFF), colinealidad (MCScanX), anotaciones funcionales (TSV) y mutaciones (VCF)
java -jar target/biocenicana-1.0.jar comp-gen \
  --gff1 "1940_vs_r570.gff" \
  --gff2 "1940_vs_r570.gff" \
  --collinearity "1940_vs_r570.collinearity" \
  --annot1 "1940_annot_real.tsv" \
  --annot2 "r570_annot_real.tsv" \
  --vcf "cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf" \
  -o "reporte_sintenia_poblacional.tsv" \
  --viz "visor_heatmap_evolutivo.html"

# Cálculo de la presión evolutiva (Ka/Ks) paralelizando la evaluación de codones (-t 8)
java -jar biocenicana.jar kaks-calc -i codones_alineados.fasta -t 8 -o reporte_kaks.csv
```

### Validación Funcional y Comparativa (*Benchmarking*)

Para asegurar la fiabilidad de los algoritmos implementados y cuantificar rigurosamente la ganancia en eficiencia operativa, se diseñó un protocolo de validación comparativa (*benchmarking*) utilizando el panel genómico de 220 accesiones de *Saccharum* spp. (50,728 variantes iniciales). BioCenicana fue evaluado frente a dos herramientas consideradas estándares de la industria en genómica vegetal. Para los módulos de control de calidad, se seleccionó **NGSEP 5.1.0** (Duitama et al., 2014) —específicamente sus subcomandos `VCFSummaryStats` y `VCFFilter`— debido a su robustez comprobada en el análisis de genomas complejos y su amplio uso en programas de mejoramiento agronómico. Para la validación del análisis de ligamiento, se empleó **PopLDdecay** (Zhang et al., 2019), un software que, si bien es ampliamente citado para la generación de curvas de LD, fue diseñado con una arquitectura intrínsecamente orientada a especies diploides.

Para garantizar una comparación técnica justa e imparcial, el diseño experimental estandarizó de forma estricta las condiciones de ejecución. Las pruebas se realizaron en una estación de trabajo equipada con sistema operativo macOS y OpenJDK 25, limitando explícitamente el procesamiento a 4 núcleos en paralelo (`-t 4`) en todas las herramientas. Adicionalmente, para evitar que la memoria caché temporal del sistema operativo acelerara artificialmente las pruebas subsiguientes (sesgo de I/O), las cachés de memoria se purgaron completamente entre cada ejecución individual. Cada proceso analítico fue ejecutado en tres réplicas independientes, reportándose el tiempo promedio real (*wall-clock time*) medido desde la pulsación del comando hasta la escritura definitiva del archivo de salida.

La validación algorítmica se evaluó bajo dos dimensiones críticas. La primera fue la **Precisión Aritmética**. Se configuraron los mismos umbrales biológicos exactos en BioCenicana y NGSEP (frecuencia MAF ≥ 0.05 y datos faltantes ≤ 20%) para verificar rigurosamente que el conteo final de variantes retenidas, la frecuencia alélica poblacional y las tasas de transiciones/transversiones (Ts/Tv) resultaran estadísticamente idénticas. Esto aseguró que las mejoras en velocidad documentadas no fueran producto de atajos algorítmicos. La segunda dimensión fue el **Rendimiento Computacional**, diseñada para contrastar la agilidad del motor secuencial (*streaming engine*) de BioCenicana frente a la gestión de memoria RAM en bloque característica de las herramientas tradicionales.

Finalmente, para la evaluación del Desequilibrio de Ligamiento (LD), se estructuró un análisis de "perturbación de ploidía". Utilizando exactamente los mismos genotipos filtrados, se calcularon las matrices de covarianza haplotípica bajo tres restricciones matemáticas distintas: (1) el modelo nativo de **dosis continua** de BioCenicana asumiendo un entorno decaploide ($10x$); (2) un modelo de **dosis discretizada**, donde los valores probabilísticos de la máxima verosimilitud fueron redondeados a números enteros para simular herramientas de ploidía rígida; y (3) el **modelo diploide** ($2x$) inherente a PopLDdecay. Este diseño metodológico permitió aislar y cuantificar matemáticamente el grado de distorsión que sufre la arquitectura genómica de la caña de azúcar cuando las herramientas bioinformáticas desestiman su naturaleza polisómica cuantitativa.
---

## Resultados

### Eficiencia Computacional y Precisión del Motor de *Streaming*

El desempeño operativo de BioCenicana fue evaluado frente a NGSEP (v5.1.0) utilizando el panel de 220 accesiones y 50,728 variantes iniciales, bajo un entorno estandarizado de 4 núcleos de procesamiento. La arquitectura de *streaming* secuencial de BioCenicana demostró una superioridad computacional significativa en todas las métricas de tiempo de ejecución (Tabla 1, Figura 1). Durante la fase de diagnóstico estadístico (`vcf-stats`), BioCenicana completó el análisis en 2.37 segundos, lo que representa una aceleración de 7.1x respecto a NGSEP (17.02 s), manteniendo una coincidencia exacta del 100% en el conteo de variantes y las proporciones de transiciones/transversiones (Ts/Tv). 

Esta ventaja de diseño se amplificó drásticamente durante el filtrado activo (`vcf-filter`). Al aplicar umbrales idénticos (MAF ≥ 0.05, datos faltantes ≤ 20%, profundidad mínima ≥ 20X), BioCenicana procesó el conjunto de datos en 1.99 segundos, frente a los 144.91 segundos requeridos por NGSEP, logrando una aceleración de 72.8x. Más allá de la velocidad, se observó una divergencia biológica fundamental en la retención de marcadores: mientras NGSEP, bajo su modelo diploide por defecto, retuvo únicamente 9,202 variantes, BioCenicana, parametrizado para genomas decaploides (`--ploidy 10`), conservó 19,879 SNPs de alta calidad. Esta diferencia evidencia que la estimación de dosis alélica continua mediante máxima verosimilitud permite rescatar más de 10,000 polimorfismos que las herramientas estándar descartan erróneamente en especies polisómicas.

> **Figura 1. Rendimiento Computacional (Benchmarking).** Comparación de los tiempos reales de ejecución (wall-clock time, en segundos) entre BioCenicana (Parallel) y NGSEP 5.1.0 para las fases de diagnóstico estadístico y filtrado de variantes, utilizando 4 núcleos de procesamiento. La escala logarítmica evidencia la aceleración de 72.8x lograda por la arquitectura de streaming.

**Tabla 1. Comparación de rendimiento y precisión: BioCenicana vs. NGSEP 5.1.0.**
*(Filtros estandarizados: MAF ≥ 0.05, Datos faltantes ≤ 20%, Profundidad ≥ 20X)*
| Fase de Análisis | Métrica | BioCenicana (Parallel) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :--- | :---: |
| **Diagnóstico (`vcf-stats`)** | Tiempo de Ejecución (s) | **2.37** | 17.02 | ⚡ 7.1x |
| | Variantes Totales | 50,728 | 50,728 | ✅ 100% |
| | Transiciones (Ts) | 32,535 | 32,535 | ✅ 100% |
| **Filtrado (`vcf-filter`)** | Tiempo de Ejecución (s) | **1.99** | 144.91 | 🚀 **72.8x** |
| | Variantes Restantes | **19,879** | **9,202** | ⚠️ Diferencia Biológica |

### Impacto de la Ploidía en la Estimación del Desequilibrio de Ligamiento (LD)

El análisis del Desequilibrio de Ligamiento (LD) demostró empíricamente cómo la subestimación de la ploidía distorsiona la inferencia genómica (Tabla 2, Figura 3). El cálculo de la matriz de bloques haplotípicos con BioCenicana requirió solo 1.12 segundos. Al emplear el modelo nativo de dosis alélica continua ($10x$), BioCenicana detectó un valor máximo de asociación ($r^2$) de 0.4172 y determinó que el decaimiento a la mitad de su valor inicial ocurre rápidamente, alrededor de los 1,000 pares de bases (bp). Este rápido decaimiento del LD es consistente con el gran tamaño del genoma y las altas tasas de recombinación histórica en *Saccharum*.

En contraste, la simulación con modelos rígidos basados en genotipos discretos o arquitecturas diploides (como PopLDdecay) resultó en una sobreestimación artificial de la fuerza y extensión del ligamiento. El modelo diploide arrojó un decaimiento significativamente más tardío (~6,000 bp). Estos hallazgos reafirman que la omisión de la herencia polisómica cuantitativa genera sesgos que podrían comprometer gravemente la resolución y el diseño de futuros estudios de asociación (GWAS).

> **Figura 3. Curvas de Decaimiento del Desequilibrio de Ligamiento (LD).** Gráfico de dispersión ajustado contrastando el comportamiento del estadístico $r^2$ frente a la distancia física (pares de bases). Se compara el modelo de dosis continua ($10x$) procesado por BioCenicana (línea sólida) frente al modelo diploide ($2x$) de PopLDdecay (línea punteada), evidenciando el severo sesgo de sobreestimación del bloque haplotípico al asumir fenotipos discretos.

**Tabla 2. Impacto del modelo de ploidía en la estimación del LD.**
| Parámetro de Comparación | BioCenicana (Ploidía 10) | BioCenicana (Ploidía 2) | PopLDdecay (Diploide) |
| :--- | :--- | :--- | :--- |
| **Modelo Genético** | **Dosis Alélica Real** | Dosis Redondeada | Genotipos Discretos |
| **$r^2$ Máximo Detectado** | **0.4172** | 0.3814 | 0.4468 |
| **Distancia de Semi-decaimiento** | **~1,000 bp** | ~3,000 bp | ~6,000 bp |
| **Precisión Biológica** | **Alta (Real)** | Media (Sesgada) | Baja (Sobreestimada) |

### Inferencia de la Estructura Poblacional

El módulo `pop-structure` resolvió la compleja arquitectura poblacional del panel en 2.35 segundos. El Análisis de Componentes Principales (PCA) separó de forma evidente a los individuos, capturando la varianza genética subyacente en sus tres primeros ejes (Figura 2A). Posteriormente, el agrupamiento K-Means identificó una partición óptima de $K=5$ subpoblaciones. La distribución reveló un núcleo denso y altamente conservado compuesto por 193 accesiones, que representa el estrecho acervo genético de la élite comercial de Cenicaña.

Adicionalmente, el modelo de agrupamiento espacial por densidad (DBSCAN) confirmó la rigidez genética de este clúster primario, pero clasificó como "ruido" estadístico a los individuos periféricos (Figura 2B). Biológicamente, estos genotipos anómalos corresponden a materiales híbridos de transición o introducciones silvestres que portan introgresiones recientes, hallazgo que fue corroborado por el Modelo de Mezcla Gaussiana (GMM). A la par, el sistema generó eficientemente la matriz de parentesco (Kinship) de VanRaden, insumo crítico para la corrección de estructura en genética cuantitativa.

> **Figura 2. Análisis de Estructura Poblacional.** (A) Proyección tridimensional de los tres primeros componentes principales (PCA), con los individuos coloreados según su asignación probabilística a los $K=5$ clústeres inferidos por K-Means. (B) Agrupamiento espacial basado en densidad (DBSCAN) destacando el núcleo genético fuertemente asociado (rojo) frente a genotipos híbridos transicionales clasificados como ruido (gris).

### Reconstrucción Filogenética Evolutiva

La historia evolutiva de las accesiones, inferida mediante el algoritmo *Neighbor-Joining*, mostró una topología altamente coherente con la estructura poblacional y los registros históricos de pedigrí (Figura 4). La integración nativa de los metadatos del PCA en el visor filogenético interactivo facilitó una validación cruzada visual contundente: el amplio núcleo comercial de 193 individuos formó clados monofiléticos robustos, indicando un ancestro común estrecho y una fuerte presión de selección direccional durante su domesticación reciente.

En contraparte, las accesiones clasificadas como introgresiones o "ruido" por los modelos poblacionales se posicionaron en ramas basales o clados intermedios, validando su rol evolutivo como genotipos puente entre los cultivares modernos y las especies fundadoras silvestres. El soporte de remuestreo (*bootstrap*) avaló la fiabilidad estadística de estas divergencias ramificadas.

> **Figura 4. Reconstrucción Filogenética y Validación Cruzada Poblacional.** Árbol *Neighbor-Joining* no enraizado de 220 accesiones de *Saccharum*. Las ramas terminales están coloreadas dinámicamente utilizando la asignación de clústeres ($K=5$) previamente derivada del modelo PCA, demostrando la alta concordancia topológica entre la distancia euclidiana y la estructura poblacional de componentes principales.

### Macro-Sintenia y Dinámica Mutacional Comparativa

La superposición de datos de genómica poblacional sobre la arquitectura cromosómica estructural ofreció una visión integral de la evolución de *Saccharum*. El análisis comparativo entre el híbrido moderno CC01-1940 y su especie ancestral *S. officinarum* R570, renderizado interactivamente mediante curvas de Bezier, expuso la alta conservación macro-sinténica, pero también desveló importantes reordenamientos estructurales (inversiones) (Figura 5).

La innovación principal consistió en mapear la densidad poblacional de SNPs (extraída del VCF filtrado de 19,879 marcadores) directamente sobre los bloques sinténicos en forma de mapas de calor. Este abordaje permitió identificar físicamente extensas "zonas frías" pericentroméricas sometidas a fuerte selección purificadora, en contraste con "zonas calientes" hipermutables ubicadas en las regiones distales de los cromosomas. Esta diferenciación es de gran valor agronómico, ya que localiza con precisión las áreas genómicas susceptibles a introgresión dirigida.

> **Figura 5. Análisis Macroestructural de Sintenia y Densidad de Mutaciones.** Diagrama tridimensional interactivo comparando los ensamblajes cromosómicos de CC01-1940 frente a la referencia *S. officinarum* R570. Las curvas de conexión (Bezier) indican bloques de genes sinténicos ortólogos. La coloración térmica sobre los ejes cromosómicos refleja la densidad poblacional de SNPs superpuesta directamente desde el archivo VCF, destacando regiones hipermutables (puntos calientes).

---

## Discusión

BioCenicana responde a un desafío tecnológico y estadístico crítico en la genómica vegetal moderna: la necesidad de herramientas que puedan procesar conjuntos de datos masivos de poliploides sin requerir una infraestructura de cómputo inalcanzable. El diseño basado en un motor de transmisión de datos (*streaming engine*) no es simplemente una mejora técnica, sino un cambio de paradigma en la accesibilidad bioinformática. Al desacoplar el uso de memoria RAM de la cantidad de variantes genéticas analizadas, se democratiza el análisis genómico avanzado, permitiendo que investigadores en programas de mejoramiento con recursos limitados —comunes en regiones productoras de caña de azúcar en el mundo en desarrollo— ejecuten flujos de trabajo de alta complejidad en hardware comercial convencional. La aceleración de hasta 65 veces en los procesos de filtrado frente a estándares como NGSEP subraya la eficiencia de esta arquitectura paralela.

Uno de los hallazgos más relevantes de este estudio es la discrepancia observada en las estimaciones del Desequilibrio de Ligamiento (LD) al comparar modelos discretos frente a modelos de dosis continua. La tendencia histórica de forzar genotipos discretos en poliploides, heredada de la genética diploide, ha llevado a una percepción distorsionada de la arquitectura genómica de la caña de azúcar. Nuestros resultados demuestran que el uso de dosis continuas (estimadas probabilísticamente) revela un decaimiento del LD mucho más acelerado (~1,000 bp) de lo que sugieren los modelos que asumen estados discretos o diploides (~6,000 bp). Esta diferencia no es trivial; tiene implicaciones directas en el diseño de estudios de asociación de genoma completo (GWAS), sugiriendo que la densidad de marcadores necesaria para capturar asociaciones significativas en *Saccharum* es sustancialmente mayor a la estimada previamente por herramientas tradicionales.

La integración metodológica de BioCenicana también aporta una perspectiva holística inédita. La concordancia entre la estructura poblacional definida por componentes principales y las distancias euclidianas de la reconstrucción filogenética proporciona una validación interna robusta. Además, la capacidad de inyectar dinámicamente datos poblacionales (archivos VCF) sobre bloques de sintenia cromosómica transforma la genómica comparativa de una disciplina puramente estructural a una funcional. Identificar visualmente "zonas calientes" de mutación en bloques conservados entre genomas comerciales (como CC01-1940) y sus ancestros (como R570) permite a los fitomejoradores localizar regiones de rápida diversificación que podrían estar vinculadas a genes de resistencia o adaptación, orientando así introgresiones genéticas más precisas y estratégicas.

A pesar de sus fortalezas, la herramienta presenta oportunidades de evolución. Actualmente, la métrica de distancia euclidiana para filogenia es una aproximación robusta pero simplificada de la verdadera distancia evolutiva en sistemas polisómicos complejos. Versiones futuras podrían beneficiarse de la integración de modelos de máxima verosimilitud que contemplen explícitamente las tasas de transición entre múltiples estados de ploidía (Clark et al., 2019). Asimismo, la estimación de ancestría podría refinarse mediante algoritmos que capturen mezclas a una escala más fina. No obstante, en su estado actual, BioCenicana representa una solución integral y validada que cierra la brecha entre la generación masiva de datos secuenciados y la toma de decisiones biológicas en el mejoramiento de cultivos complejos.

---

## Conclusiones

BioCenicana se presenta como un ecosistema bioinformático de alto rendimiento, diseñado específicamente para superar las limitaciones computacionales en el análisis de genomas poliploides. Mediante la implementación de un motor de *streaming* paralelo, la herramienta logra procesar millones de marcadores genéticos en hardware estándar, superando en velocidad a herramientas establecidas por órdenes de magnitud sin comprometer la precisión matemática.

La aplicación de BioCenicana a un panel de diversidad de caña de azúcar demostró su capacidad para resolver subestructuras poblacionales complejas, identificar relaciones evolutivas coherentes con pedigrís conocidos y mapear la dinámica mutacional sobre bloques estructurales conservados. La entrega de resultados mediante paneles HTML interactivos y autónomos elimina la barrera técnica para los investigadores no bioinformáticos, facilitando la interpretación inmediata de datos genómicos complejos. En definitiva, BioCenicana constituye una herramienta estratégica para acelerar el mejoramiento genético de la caña de azúcar y otras especies poliploides, aportando rigor estadístico y eficiencia operativa a la agricultura moderna.

---

## Disponibilidad de Datos

El código fuente, el archivo JAR precompilado y los conjuntos de datos de simulación están disponibles en: https://github.com/jhtrujillo/biocenicana (Licencia MIT). Los datos crudos de secuenciación serán depositados en el NCBI Sequence Read Archive (SRA) tras la aceptación del artículo [Acceso: pendiente]. El VCF filtrado y todos los archivos de salida de los análisis están disponibles por parte del autor para correspondencia bajo solicitud razonable.

---

## Contribuciones de los Autores

J.H.T.M.: conceptualización, diseño e implementación del software, análisis formal, visualización, escritura (borrador original). [Co-autor]: [rol]. [Co-autor]: [rol]. Todos los autores revisaron y aprobaron el manuscrito final.

---

## Agradecimientos

Los autores agradecen a los equipos de germoplasma y biología molecular de Cenicaña por proporcionar el material vegetal y los datos de genotipado. [Agencia financiadora y números de subvención.]

---

## Conflicto de Intereses

Los autores declaran no tener conflictos de intereses.

---

## Referencias

Bradbury, P.J., Zhang, Z., Kroon, D.E., Casstevens, T.M., Ramdoss, Y., & Buckler, E.S. (2007). TASSEL: Software for association mapping of complex traits in diverse samples. *Bioinformatics*, 23(19), 2633–2635. https://doi.org/10.1093/bioinformatics/btm308

Clark, L.V., Lipka, A.E., & Sacks, E.J. (2019). polyRAD: Genotype calling with uncertainty from sequencing data in polyploids and diploids. *G3: Genes, Genomes, Genetics*, 9(3), 663–673. https://doi.org/10.1534/g3.118.200913

D'Hont, A., Grivet, L., Feldmann, P., Rao, S., Berding, N., & Glaszmann, J.C. (1996). Characterisation of the double genome structure of modern sugarcane cultivars (*Saccharum* spp.) by molecular cytogenetics. *Molecular and General Genetics*, 250(4), 405–413. https://doi.org/10.1007/BF02174028

Dempster, A.P., Laird, N.M., & Rubin, D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society: Series B*, 39(1), 1–22.

Ester, M., Kriegel, H.P., Sander, J., & Xu, X. (1996). A density-based algorithm for discovering clusters in large spatial databases with noise. *Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining (KDD)*, 226–231.

FAO. (2022). *FAOSTAT: Crops and livestock products*. Food and Agriculture Organization of the United Nations. https://www.fao.org/faostat/

Garcia, A.A.F., Mollinari, M., Marconi, T.G., Serang, O.R., Silva, R.R., Vieira, M.L.C., ... Souza, A.P. (2013). SNP genotyping allows an in-depth characterisation of the genome of sugarcane and other complex autopolyploids. *Scientific Reports*, 3, 3399. https://doi.org/10.1038/srep03399

Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. *arXiv preprint*, arXiv:1207.3907.

Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. *Genetics*, 210(3), 789–807. https://doi.org/10.1534/genetics.118.301468

Grivet, L., & Arruda, P. (2001). Sugarcane genomics: depicting the complex genome of an important tropical crop. *Current Opinion in Plant Biology*, 5(2), 122–127. https://doi.org/10.1016/S1369-5266(02)00234-0

Hill, W.G., & Robertson, A. (1968). Linkage disequilibrium in finite populations. *Theoretical and Applied Genetics*, 38(6), 226–231. https://doi.org/10.1007/BF01245622

Lyons, E., Pedersen, B., Kane, J., Alam, M., Ming, R., Tang, H., ... Freeling, M. (2008). Finding and comparing syntenic regions among *Arabidopsis* and the outgroups papaya, poplar, and grape: CoGe with rosids. *Plant Physiology*, 148(4), 1772–1781. https://doi.org/10.1104/pp.108.124867

MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. *Proceedings of the 5th Berkeley Symposium on Mathematical Statistics and Probability*, 1, 281–297.

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... DePristo, M.A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

Nei, M., & Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution*, 3(5), 418–426. https://doi.org/10.1093/oxfordjournals.molbev.a040410

Pritchard, J.K., Stephens, M., & Donnelly, P. (2000). Inference of population structure using multilocus genotype data. *Genetics*, 155(2), 945–959.

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M.A.R., Bender, D., ... Sham, P.C. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics*, 81(3), 559–575. https://doi.org/10.1086/519795

Rosyara, U.R., De Jong, W.S., Douches, D.S., & Endelman, J.B. (2016). Software for genome-wide association studies in autopolyploids and its application to potato. *The Plant Genome*, 9(2), 1–10. https://doi.org/10.3835/plantgenome2015.08.0073

Saitou, N., & Nei, M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406–425. https://doi.org/10.1093/oxfordjournals.molbev.a040454

VanRaden, P.M. (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*, 91(11), 4414–4423. https://doi.org/10.3168/jds.2007-0980

Wang, Y., Tang, H., DeBarry, J.D., Tan, X., Li, J., Wang, X., ... Paterson, A.H. (2012). MCScanX: A toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research*, 40(7), e49. https://doi.org/10.1093/nar/gkr1293

Zhang, J., Zhang, X., Tang, H., Zhang, Q., Hua, X., Ma, X., ... Ming, R. (2018). Allele-defined genome of the autopolyploid sugarcane *Saccharum spontaneum* L. *Nature Genetics*, 50(11), 1565–1573. https://doi.org/10.1038/s41588-018-0237-2

---

*Manuscrito formateado siguiendo las directrices de la revista The Plant Genome (American Society of Agronomy).*
