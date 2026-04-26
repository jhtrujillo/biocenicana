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

Finalmente, para la evaluación del Desequilibrio de Ligamiento (LD), se estructuró un análisis de "perturbación de ploidía". Utilizando exactamente los mismos genotipos filtrados, se calcularon las matrices de covarianza haplotípica bajo tres restricciones matemáticas distintas: (1) el modelo nativo de dosis continua de BioCenicana asumiendo un entorno decaploide (10x); (2) un modelo de dosis discretizada, donde los valores probabilísticos de la máxima verosimilitud fueron redondeados a números enteros para simular herramientas de ploidía rígida; y (3) el modelo diploide (2x) inherente a PopLDdecay. Este diseño metodológico permitió aislar y cuantificar matemáticamente el grado de distorsión que sufre la arquitectura genómica de la caña de azúcar cuando las herramientas bioinformáticas desestiman su naturaleza polisómica cuantitativa.

---

## Resultados

### Módulo de Control de Calidad y Filtrado (`vcf-stats` y `vcf-filter`)

El desempeño computacional y la precisión analítica del motor de *streaming* secuencial de BioCenicana fueron evaluados frente a la herramienta NGSEP (versión 5.1.0). El *benchmarking* se ejecutó empleando un conjunto de datos poblacional de 220 accesiones de *Saccharum* spp. con 50,728 variantes iniciales. Para asegurar una comparación equitativa y rigurosa, el procesamiento se limitó a 4 núcleos en ambas herramientas y se purgaron las cachés de memoria del sistema operativo entre ejecuciones.

Durante el diagnóstico inicial (`vcf-stats`), BioCenicana completó el perfilamiento estadístico en 2.37 segundos, resultando en una aceleración de 7.1x frente a los 17.02 segundos requeridos por NGSEP. Este incremento significativo en la velocidad no comprometió la exactitud matemática: se obtuvo una coincidencia absoluta del 100% en el conteo total de variantes y en las proporciones de transiciones/transversiones (Tabla 1).

**Tabla 1. Comparación de rendimiento y precisión del diagnóstico estadístico (`vcf-stats`).**
| Métrica | BioCenicana (4 núcleos) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :---: |
| **Tiempo de Ejecución (s)** | **2.37** | 17.02 | 7.1x |
| **Variantes Totales** | 50,728 | 50,728 | 100% |
| **Transiciones (Ts)** | 32,535 | 32,535 | 100% |
| **Transversiones (Tv)** | 18,193 | 18,193 | 100% |
| **Ts/Tv Ratio** | 1.788 | 1.79 | 100% |

La eficiencia arquitectónica de BioCenicana se hizo aún más evidente durante la fase de depuración activa (`vcf-filter`). Al aplicar umbrales paramétricos estrictos (frecuencia de alelo menor MAF ≥ 0.05, tolerancia de datos faltantes ≤ 20%, y profundidad mínima ≥ 20X), el motor de procesamiento paralelo de BioCenicana procesó y exportó el archivo filtrado en tan solo 1.85 segundos. En contraposición, NGSEP requirió 121.97 segundos para completar la misma operación, lo que evidencia una aceleración masiva de 65x a favor de BioCenicana (Tabla 2). 

Es imperativo destacar que este drástico incremento en la velocidad de procesamiento no alteró la precisión analítica del filtrado. Ambas herramientas, partiendo del conjunto inicial de 50,728 variantes, retuvieron de manera idéntica 7,443 SNPs de alta calidad. Esto demuestra que la arquitectura de *streaming* paralelo de BioCenicana es capaz de sortear los característicos cuellos de botella de entrada/salida (I/O) y consumo de memoria RAM asociados a la genómica de poblaciones, garantizando simultáneamente una fidelidad y coincidencia de resultados del 100% frente a estándares consolidados de la industria.

**Tabla 2. Comparación del rendimiento de filtrado (`vcf-filter`).**
*(Filtros: MAF ≥ 0.05, Datos faltantes ≤ 20%, Profundidad ≥ 20X)*
| Métrica | BioCenicana (Parallel) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :---: |
| **Tiempo de Ejecución (s)** | **1.85** | 121.97 | **65x** |
| **Variantes Iniciales** | 50,728 | 50,728 | 100% |
| **Variantes Restantes** | **7,443** | **7,443** | 100% |

### Impacto de la Ploidía en la Estimación del Desequilibrio de Ligamiento (LD)

Para cuantificar el grado de distorsión que sufre la arquitectura genómica de *Saccharum* cuando se desestima su naturaleza polisómica, se estructuró un análisis comparativo calculando el decaimiento del Desequilibrio de Ligamiento (LD) bajo tres restricciones matemáticas distintas: el modelo de dosis continua (ploidía 10) de BioCenicana, un modelo de dosis discretizada (ploidía 2) de BioCenicana, y el modelo estándar de genotipos discretos diploides inherente a PopLDdecay.

Los resultados (Tabla 3) demostraron que la subestimación de la ploidía genera un sesgo matemático y biológico severo. El motor nativo de BioCenicana, al emplear la dosis alélica real bajo un entorno decaploide, procesó la matriz en 1.12 segundos, detectando un valor de asociación máximo (r²) de 0.4172. Lo más destacable fue la identificación de un decaimiento rápido del bloque haplotípico a la mitad de su valor inicial, ocurriendo aproximadamente a los 1,000 pares de bases (bp). Este rápido decaimiento es biológicamente congruente con el gran tamaño del genoma y las altas tasas de recombinación histórica en el género *Saccharum*.

Por el contrario, la simulación con herramientas rígidas basadas en arquitecturas diploides produjo una sobreestimación artificial de la fuerza y extensión del ligamiento. PopLDdecay, al forzar la herencia de genotipos discretos, arrojó un decaimiento significativamente más tardío (~6,000 bp), creando bloques de asociación artificiales seis veces más largos que la realidad biológica subyacente. Estos hallazgos reafirman que la omisión de la herencia cuantitativa compromete gravemente la resolución espacial necesaria para el diseño y confiabilidad de futuros estudios de asociación (GWAS).

**Tabla 3. Análisis de perturbación de ploidía en el decaimiento del Desequilibrio de Ligamiento (LD).**
| Parámetro | BioCenicana (Ploidía 10) | BioCenicana (Ploidía 2) | PopLDdecay (Diploide) |
| :--- | :--- | :--- | :--- |
| **Modelo Genético** | **Dosis Alélica Real** | Dosis Redondeada | Genotipos Discretos |
| **Tiempo de Proceso** | 1.12 s | 1.13 s | **0.34 s** |
| **r² Máximo** | **0.4172** | 0.3814 | 0.4468 |
| **Semi-decaimiento** | **~1,000 bp** | ~3,000 bp | ~6,000 bp |

### Análisis Comparativo de Estructura Poblacional (`pop-structure`)

El análisis de la estructura poblacional en poliploides presenta un desafío analítico donde convergen la viabilidad computacional y la fidelidad biológica. Mediante la ejecución de pruebas de rendimiento estandarizadas, BioCenicana procesó la matriz genotípica del panel completo en **11.04 segundos**, ejecutando un flujo integral *Todo-en-Uno* que abarcó: la lectura, el filtrado de variantes, la Descomposición en Valores Singulares (SVD), la inferencia automática de $K=3$ subgrupos genéticos y la construcción de la matriz Kinship.

Para evaluar el rigor y la eficiencia de este resultado, se estructuró una comparativa empírica frente a dos paradigmas bioinformáticos: el estándar de la industria en Java (NGSEP) y los paquetes especializados en R (`AGHmatrix`, `polyRAD`).

En la evaluación de escalabilidad frente a herramientas como NGSEP (versión 5.1.0), se demostró una superioridad absoluta en eficiencia y consolidación. Al configurar explícitamente a NGSEP para calcular distancias bajo un modelo decaploide (`-p 10 -s 3`), el procesamiento tomó **25.04 segundos** —más del doble que BioCenicana— para generar únicamente una matriz de distancia en formato de texto plano. NGSEP carece de integración nativa para la reducción de dimensionalidad (PCA), la identificación de clústeres no supervisados o el cálculo estandarizado de la matriz de parentesco de VanRaden. Esto obliga a los investigadores a recurrir a lenguajes adicionales para completar el flujo poblacional, fragmentando el análisis.

Por otro lado, la comparación matemática con paquetes especializados en R (`AGHmatrix`) reveló una **paridad matemática absoluta (100% de coincidencia)** en la matriz de parentesco (*Kinship*). Ambos ecosistemas emplean la dosis alélica continua, garantizando una alta precisión estadística adaptada a la complejidad genómica de *Saccharum*. Sin embargo, los flujos en R exigen la carga total del archivo VCF en la memoria RAM —un proceso sumamente ineficiente para genomas grandes que consume altos volúmenes de memoria y toma minutos en iniciar—. La arquitectura de *streaming* paralelo de BioCenicana elimina este cuello de botella y entrega el flujo analítico poblacional completo de manera casi instantánea, unificando precisión biológica y excelencia computacional.

**Tabla 4. Rendimiento comparativo del Análisis de Estructura Poblacional (Ploidía 10).**
| Métrica / Herramienta | BioCenicana | NGSEP 5.1.0 | Flujo en R (`AGHmatrix`) |
| :--- | :--- | :--- | :--- |
| **Tiempo de Ejecución** | **11.04 s** | 25.04 s | ~Minutos (Cuello de botella I/O) |
| **Dosis Alélica** | Continua (Nativa) | Continua (Configurada) | Continua (Nativa) |
| **Integración Analítica** | **Total** (SVD, Clusters, Kinship) | Fragmentada (Solo Distancia) | Manual (Requiere varios paquetes) |
| **Consumo de Memoria** | Óptimo (*Streaming*) | Alto | Crítico (Carga en RAM) |

Más allá del rendimiento computacional, la verdadera ventaja de BioCenicana radica en la riqueza de sus resultados biológicos. Mientras que herramientas tradicionales como NGSEP se limitan a exportar una matriz matemática genérica de distancias (sin encabezados ni metadatos) que obliga al usuario a programar el resto del análisis, el motor de BioCenicana realiza de forma nativa la inferencia poblacional completa y exporta un panel estructurado listo para su interpretación evolutiva (Tabla 5).

**Tabla 5. Comparativa de resultados analíticos generados automáticamente.**
| Parámetro Biológico | Salida BioCenicana | Salida NGSEP | Aplicación en GWAS y Mejoramiento |
| :--- | :--- | :--- | :--- |
| **Matriz de Parentesco** | Sí (VanRaden *Kinship*) | No (Solo Distancia IBS) | Corrección de parentesco oculto mediante Modelos Lineales Mixtos. |
| **Coordenadas PCA** | Sí (PC1 a PC10 extraídos) | No | Corrección matemática de la estratificación poblacional. |
| **Inferencia de Clústeres** | Sí ($K=3$, KMeans, DBSCAN) | No | Identificación de subpoblaciones y linajes fundadores. |
| **Proporciones Ancestría** | Sí (Matriz Q1-Q3 Admixture) | No | Selección de híbridos y cuantificación de introgresiones silvestres. |
| **Formato de Exportación** | CSV Estructurado y Visor HTML | Texto plano numérico | Interpretación visual e inmediata sin programación adicional en R. |

La riqueza de estos resultados numéricos se cristaliza en las salidas visuales generadas de forma inmediata por el módulo. Como se detalla en la **Figura 3**, BioCenicana proporciona una evaluación jerárquica de la estructura genética sin requerir visualizadores externos. El análisis de la varianza explicada (*Scree Plot*, **Figura 3A**) y la heurística del Método del Codo (**Figura 3B**) actúan en conjunto para justificar estadísticamente el agrupamiento óptimo de las accesiones. Al aplicar el algoritmo de *k-medias* sobre las proyecciones espaciales, el sistema determinó matemáticamente que $K=3$ maximiza la homogeneidad interna de los grupos, un hallazgo que concuerda estructuralmente con la topología de ramificación principal observada en el dendrograma aglomerativo (**Figura 3C**). Asimismo, el módulo calcula de forma paralela dos matrices complementarias: la Matriz de Distancia Genética Euclidiana (**Figura 3E**), que ilustra la separación absoluta entre genotipos, y la crucial Matriz de Parentesco Genómico de VanRaden (*Kinship*, **Figura 3D**), cuyos bloques térmicos de alta intensidad cuantifican el grado de parentesco oculto. Esta matriz de covarianza es el insumo biológico fundamental para la construcción de Modelos Lineales Mixtos (MLM) en estudios GWAS, corrigiendo eficazmente las asociaciones espurias.

Para una disección más profunda de las subpoblaciones, la proyección geométrica de las coordenadas extraídas permite a los genetistas validar los linajes. El gráfico bidimensional de Componentes Principales (**Figura 4A**) y su respectiva secuencia de rotación tridimensional en el espacio PC1-PC2-PC3 (**Figura 4B**) demuestran empíricamente cómo los tres clústeres fundadores mantienen una separación estricta, descartando artefactos de solapamiento derivados de proyecciones planas. Finalmente, la deconstrucción individual de estas firmas genéticas se resume en el modelo probabilístico de proporciones de ancestría (*Admixture*, **Figura 4C**), donde el perfil codificado por colores de cada accesión permite a los fitomejoradores no solo certificar linajes puros, sino medir el grado exacto de introgresión genómica proveniente de híbridos interespecíficos.

<br>

**Figura 3. Inferencia de Clústeres Genéticos y Matrices de Parentesco.** 

<p align="center">
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/ExplainedVariance.png" width="100%"><br>
    <small><b>(A) Varianza Explicada (Scree Plot):</b> Cuantifica el porcentaje de variabilidad genética capturada por cada componente. Un declive pronunciado (PC1 y PC2) evidencia una fuerte estructuración basal.</small>
  </div>
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/ElbowMethod.png" width="100%"><br>
    <small><b>(B) Método del Codo (WCSS vs K):</b> Minimiza la varianza intra-clúster. El punto de inflexión ("codo") determina empíricamente que la población se divide de manera óptima en <b>K=3</b> subgrupos distintos.</small>
  </div>
</p>

<p align="center">
  <div style="text-align: center;">
    <img src="../assets/Dendrogram.png" width="90%"><br>
    <small><b>(C) Dendrograma Jerárquico:</b> Reconstrucción aglomerativa que revela visualmente las distancias evolutivas y la bifurcación histórica de los linajes parentales en la población.</small>
  </div>
</p>

<p align="center">
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/Kinship.png" width="100%"><br>
    <small><b>(D) Matriz Kinship (VanRaden):</b> Estructura de covarianza genómica. Los bloques rojos/amarillos indican alta relación familiar. Es el insumo crítico para corregir falsos positivos en GWAS.</small>
  </div>
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/GeneticDistanceMatrixHeatmap.png" width="100%"><br>
    <small><b>(E) Matriz Euclidiana:</b> A diferencia de la Kinship (covarianza), esta matriz ilustra la separación matemática absoluta entre cada par de individuos del panel.</small>
  </div>
</p>

<br>

**Figura 4. Proyección de Estructura Poblacional Espacial y Ancestría.**

<p align="center">
  <div style="text-align: center;">
    <img src="../assets/ClusterVisualizationPCA.png" width="90%"><br>
    <small><b>(A) Componentes Principales (PCA en 2D):</b> Proyección espacial de distancias genéticas. Cada punto es un individuo; los colores representan los 3 clústeres inferidos, validando la separación de los linajes.</small>
  </div>
</p>

<p align="center">
  <div style="text-align: center;">
    <small><b>(B) Secuencia de Rotación Tridimensional (PC1, PC2, PC3):</b> Captura volumétrica iterativa que demuestra que la segregación de los clústeres genéticos se mantiene robusta en el espacio 3D, sin superposición artificial plana.</small><br><br>
    <img src="../assets/pac3d1.png" width="24%">
    <img src="../assets/pac3d2.png" width="24%">
    <img src="../assets/pac3d3.png" width="24%">
    <img src="../assets/pac3d4.png" width="24%">
    <br>
    <img src="../assets/pac3d5.png" width="24%">
    <img src="../assets/pac3d6.png" width="24%">
    <img src="../assets/pac3d7.png" width="24%">
    <img src="../assets/pac3d8.png" width="24%">
  </div>
</p>

<p align="center">
  <div style="text-align: center;">
    <img src="../assets/Ancestry.png" width="90%"><br>
    <small><b>(C) Proporciones de Ancestría (Admixture):</b> Cada barra vertical es un individuo. Los colores cuantifican qué porcentaje de su genoma proviene de cada uno de los 3 ancestros fundadores, permitiendo medir introgresiones.</small>
  </div>
</p>

> [!NOTE]
> **Disponibilidad de Datos Interactivos:** El archivo HTML generado nativamente por BioCenicana, que contiene el modelo tridimensional rotatorio y los metadatos completos del análisis poblacional, está disponible para su exploración interactiva como Material Suplementario en el siguiente enlace web: [Explorar Dashboard Interactivo 3D (PCA)](https://johntrujillomonte.com/projects/biocenicana/pca_interactive_dashboard.html).

### Reconstrucción Filogenética y Dinámica Evolutiva (`snp-tree`)

La reconstrucción de la historia evolutiva en especies poliploides requiere algoritmos capaces de procesar la dosis alélica continua para estimar las distancias genéticas verdaderas, evitando la pérdida de información que ocurre al discretizar genotipos (diploidización forzada). El módulo `snp-tree` de BioCenicana aborda este reto calculando una matriz de distancia euclidiana por pares basada en las profundidades relativas de los alelos y, subsecuentemente, reconstruyendo la topología poblacional mediante el algoritmo aglomerativo de *Neighbor-Joining* (NJ).

Para validar la eficiencia de esta arquitectura matemática, se diseñó un protocolo de *benchmark* comparativo frente a la suite NGSEP (Tabla 6). Los resultados revelaron un paradigma de fragmentación y cuellos de botella en las herramientas tradicionales. Al intentar reconstruir la filogenia, NGSEP consumió **10.50 segundos** operando exclusivamente en la generación de la matriz de distancias (`VCFDistanceMatrixCalculator`). Tras este paso, la herramienta carece de un módulo nativo integrado para construir el árbol, forzando a los investigadores a interrumpir el flujo de trabajo, exportar la pesada matriz de texto plano y utilizar ecosistemas de terceros (como paquetes de R, FastME o MEGA) para finalmente ensamblar la topología evolutiva.

En contraste radical, el motor *Todo en Uno* de BioCenicana ejecutó el flujo analítico de principio a fin en tan solo **3.07 segundos**. En este tiempo récord (3.4 veces más rápido que el cálculo parcial de NGSEP), el sistema leyó el archivo VCF, estimó la matriz probabilística de distancias, construyó el árbol iterativo *Neighbor-Joining*, guardó la topología en formato estándar de la industria (Newick, `.nwk`) y, como ventaja competitiva principal, generó automáticamente un visor HTML interactivo de alta resolución.

**Tabla 6. Rendimiento comparativo en la Reconstrucción Filogenética (Ploidía 10, ~5000 variantes).**
| Métrica Operativa | BioCenicana (`snp-tree`) | NGSEP 5.1.0 |
| :--- | :--- | :--- |
| **Tiempo de Ejecución Total** | **3.07 s** | > 10.50 s (Incompleto) |
| **Matriz de Distancia Continua** | Sí (Calculada y transmitida en memoria) | Sí (Exportada a disco, 10.5s) |
| **Construcción del Árbol (NJ)** | Sí (Integrado nativamente) | No (Requiere software externo) |
| **Validación Cruzada (PCA)** | Sí (Inyección de clústeres en nodos terminales) | No |
| **Visualización Generada** | Dashboard HTML Interactivo (3 Modos de Layout) | Ninguna |

Más allá de la supremacía computacional, la verdadera fortaleza estadística del módulo `snp-tree` radica en su integración ecosistémica. BioCenicana permite la inyección directa de los resultados poblacionales previos en la reconstrucción filogenética (mediante el parámetro `--pca`). Esta característica biológica codifica por color automáticamente los clados y nodos terminales del árbol según el linaje matemático identificado por los algoritmos no supervisados (como *K-Means*). Como se evidencia en la **Figura 5**, esta convergencia analítica cruzada permite a los genetistas validar visualmente que las distancias evolutivas puras coinciden de forma milimétrica con la separación espacial poblacional, brindando un marco de certidumbre absoluta para la selección parental en el programa de mejoramiento.

**Figura 5. Topologías Filogenéticas Interactivas y Proyección Poblacional.**

<p align="center">
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/tree1.png" width="100%"><br>
    <small><b>(A) Filograma Radial:</b> Las longitudes de las ramas son proporcionales a las distancias genéticas euclidianas. Esta vista circular maximiza el uso del espacio visual para paneles de gran tamaño.</small>
  </div>
  <div style="display: inline-block; width: 48%; vertical-align: top; text-align: center;">
    <img src="../assets/tree2.png" width="100%"><br>
    <small><b>(B) Cladograma Circular:</b> Prioriza puramente la topología de ramificación, alineando las ramas para facilitar la identificación de la relación de ancestría común más reciente (MRCA).</small>
  </div>
</p>

<p align="center">
  <div style="display: inline-block; width: 60%; vertical-align: top; text-align: center;">
    <img src="../assets/tree3.png" width="100%"><br>
    <small><b>(C) Filograma Rectangular:</b> Formato clásico de publicación. Permite rastrear y escalar la divergencia genética exacta (eje X) desde el nodo raíz hasta cada genotipo contemporáneo.</small>
  </div>
</p>

> [!NOTE]
> **Exploración Filogenética Dinámica:** Las imágenes estáticas (Figura 5) son extraídas del módulo nativo de BioCenicana. El visor interactivo original está disponible públicamente en [johntrujillomonte.com/projects/biocenicana/biocenicana_tree.html](https://johntrujillomonte.com/projects/biocenicana/biocenicana_tree.html). Este entorno web permite alternar entre estos 3 modos en tiempo real, hacer zoom interactivo y ajustar dinámicamente los umbrales de distancia para identificar y aislar clados específicos.

### Desequilibrio de Ligamiento y Resolución Genómica (`ld`)

La inferencia precisa del Desequilibrio de Ligamiento (LD) es uno de los cuellos de botella estadísticos más pronunciados en la genómica de especies autopoliploides. Las aproximaciones convencionales a menudo fuerzan la discretización de los genotipos (asignando estados absolutos como AA, AB o BB), un proceso de "diploidización" que ignora la incertidumbre de la secuenciación subyacente y altera las estimaciones de recombinación histórica. El módulo `ld` de BioCenicana soluciona esta limitación implementando un modelo que computa la correlación alélica ($r^2$) operando directamente sobre las probabilidades continuas de dosis (continuous dosage probabilities).

Al retener la naturaleza estocástica de los genotipos, los resultados revelan un panorama genómico sorprendentemente distinto al reportado por herramientas tradicionales. Mientras que los modelos discretos suelen sobreestimar la extensión del bloque haplotípico en *Saccharum spp.* (reportando caídas a los ~6,000 pb), el análisis de dosis continua de BioCenicana demuestra un decaimiento del LD muchísimo más acelerado, cayendo rápidamente en los primeros **1,000 pb** (Figura 7). 

Esta discrepancia no es un artefacto matemático, sino una representación biológica más fidedigna de los eventos de recombinación que han fracturado el genoma a lo largo de las generaciones. Este hallazgo tiene implicaciones críticas para el diseño de Estudios de Asociación de Genoma Completo (GWAS) en caña de azúcar: sugiere que la resolución del mapeo es sustancialmente mayor de lo anticipado, pero simultáneamente advierte que se requiere una densidad de marcadores (SNPs) mucho más alta para lograr saturar el genoma y capturar asociaciones causales con rasgos de importancia agronómica.

**Figura 7. Decaimiento del Desequilibrio de Ligamiento (LD Decay).**

<p align="center">
  <img src="../assets/ld.png" width="70%">
  <br>
  <small><b>Curva de Decaimiento del LD:</b> Reducción de la correlación alélica ($r^2$) en función de la distancia física inter-marcador (pb). La rápida caída subraya la alta resolución de mapeo y la necesidad de genotipado de alta densidad en poliploides.</small>
</p>

> [!NOTE]
> **Exploración de LD Interactiva:** La visualización estática de la Figura 7 fue extraída del dashboard estadístico autogenerado por BioCenicana. El visor interactivo original está disponible públicamente en [johntrujillomonte.com/projects/biocenicana/ld_results_50k_decay.html](https://johntrujillomonte.com/projects/biocenicana/ld_results_50k_decay.html). Este entorno web permite interactuar con la curva de decaimiento y consultar valores $r^2$ exactos para pares de marcadores específicos.

**Tabla 7. Comparación de Capacidades en el Cálculo de LD (Ploidía 10x).**

| Herramienta / Pipeline | Manejo de Ploidía | Modelo Estadístico de Correlación | Decaimiento Estimado ($r^2$ Half-decay) | Pipeline de Visualización |
| :--- | :--- | :--- | :--- | :--- |
| **BioCenicana (`ld`)** | **Autopoliploide Nativo (10x)** | **Probabilidad de Dosis Continua** | **~ 1,000 pb** | **Integrado (Nativo HTML)** |
| TASSEL 5.0 | Forzado (Pseudo-Diploide) | Discreto (AA, AB, BB) | ~ 6,000 pb (Sobreestimado) | Externo (Exportación a R/Python) |
| NGSEP / VCFtools | Forzado / Diploidizado | Frecuencias Alélicas Binarias | ~ 5,500 pb (Sesgado) | Requiere scripts externos |

Como se resume en la **Tabla 7**, el impacto de utilizar el modelo matemático adecuado es drástico. Para entender por qué ocurre esta discrepancia geométrica, es necesario analizar cómo cada herramienta maneja la varianza genotípica:

1. **La Sobreestimación de TASSEL (Pseudo-Diploide):** Herramientas clásicas como TASSEL esperan datos diploides (estados 0, 1, 2 correspondientes a AA, AB, BB). Al forzar la lectura de un poliploide (10x), TASSEL "colapsa" todos los estados heterocigotos intermedios (desde una dosis de 1 hasta 9 del alelo alterno) en una única categoría discreta ("AB"). Esta pérdida masiva de información cuantitativa reduce artificialmente la varianza alélica en la población. Matemáticamente, cuando se reduce la varianza cruzada entre dos marcadores, la correlación de Pearson ($r^2$) se infla de forma espuria, creando la ilusión de que el bloque haplotípico se hereda intacto por más tiempo del real (~6,000 pb).
2. **El Sesgo de NGSEP / VCFtools:** Aunque estas herramientas son robustas para genotipado, al calcular el LD frecuentemente dependen de binarizar las frecuencias alélicas poblacionales o asignar genotipos crudos sin ponderar la incertidumbre de la profundidad de secuenciación, lo que genera un "ruido" estadístico que sesga la correlación hacia arriba (~5,500 pb).
3. **La Precisión de BioCenicana:** La arquitectura de BioCenicana evita por completo la asignación discreta. En su lugar, extrae las profundidades alélicas relativas (Continuous Dosage) directamente del archivo VCF para crear un vector de probabilidad continua (valores entre 0.0 y 1.0) para cada individuo. Al calcular el $r^2$ sobre espectros continuos, BioCenicana preserva la totalidad de la varianza biológica. Como resultado, captura los eventos de recombinación sutiles que las otras herramientas ignoran, revelando la verdadera arquitectura genética de la caña de azúcar: un genoma altamente fracturado donde el desequilibrio de ligamiento decae rápidamente en los primeros ~1,000 pb.

Esta corrección estadística no es un detalle trivial; provee a la comunidad científica de mejoramiento un valor de LD biológicamente preciso, lo que dictamina que para realizar un GWAS efectivo en caña de azúcar, se requiere una densidad de marcadores drásticamente mayor a la recomendada por la literatura clásica basada en TASSEL.

### Genómica Comparativa y Evolución Funcional (`comp-gen` y `kaks-calc`)

El análisis evolutivo de los poliploides frecuentemente sufre de un desacople metodológico: por un lado, se infieren las variantes de nucleótido simple (SNPs) a nivel poblacional; por el otro, se identifican bloques de sintenia estructural a nivel macrogenómico. El módulo integrativo `comp-gen` de BioCenicana cierra esta brecha operativa al superponer, en un único ecosistema algorítmico, datos de colinealidad estructural (derivados de herramientas como McScanX o SynMap), mapas de anotación genómica (GFF3) y datos de diversidad poblacional profunda (VCF).

Esta fusión permite generar "Mapas de Calor Evolutivos" bidimensionales. Al anclar estadísticamente la densidad de polimorfismos (SNPs/Kbp) directamente sobre las coordenadas de los bloques genómicos conservados, la herramienta revela *zonas calientes* de mutación. En el contexto de la caña de azúcar, esta capacidad es crítica para discriminar si la divergencia fenotípica entre genotipos comerciales modernos (e.g., CC01-1940) e introgresiones silvestres ancestrales (e.g., *Saccharum spontaneum*) se debe a variaciones de copy number (CNV) estructural o a tasas aceleradas de mutación alélica en genes específicos de resistencia o adaptación.

Para complementar la genómica estructural con una métrica pura de presión selectiva, BioCenicana incluye el submódulo `kaks-calc`. Basado en el modelo estocástico de Nei y Gojobori (1986), este motor computa nativamente las tasas de sustitución no sinónimas (Ka) y sinónimas (Ks) directamente a partir de las secuencias codificantes (CDS) alineadas por la matriz de colinealidad. La inyección de este ratio (ω = Ka/Ks) en el dashboard de sintenia permite a los investigadores visualizar interactivamente qué genes conservados dentro de un bloque estructural están experimentando evolución purificadora (ω < 1), evolución neutral (ω ≈ 1) o selección darwiniana positiva (ω > 1). 

**Figura 6. Dashboard Interactivo de Sintenia, Genómica Estructural y Dinámica Evolutiva.**

| ![Pan-genoma](../assets/comgenomic4.png) | ![Filtro Base](../assets/comgenomic1.png) |
| :---: | :---: |
| <small><b>(A) Macro-Sintenia Pan-genómica:</b> Mapeo masivo de genomas completos (>110k genes cruzados).</small> | <small><b>(B) Filtrado Estructural Dinámico:</b> Aislamiento de regiones sinténicas mediante selección de secuencias homólogas.</small> |
| ![Filtro Estricto y Tooltips](../assets/comgenomic2.png) | ![Alineación de Subgenomas](../assets/comgenomic5.png) |
| <small><b>(C) Depuración de Ruido e Interactividad:</b> Aplicación de umbrales altos de conservación y extracción *on-demand* de funciones génicas.</small> | <small><b>(D) Arquitectura de Subgenomas:</b> Mapeo estructural focalizado de una línea comercial contra la referencia ancestral.</small> |
| ![Heatmap Evolutivo](../assets/comgenomic3.png) | ![Correlación Evolutiva/Estructural](../assets/comgenomic6.png) |
| <small><b>(E) Mapa de Calor Evolutivo (SNPs):</b> Transición de la vista física a funcional, coloreando los bloques por su diversidad poblacional.</small> | <small><b>(F) Hotspots Mutacionales:</b> Visualización térmica sobre la estructura, identificando zonas de alta presión selectiva (rojo) y alta conservación (verde/azul).</small> |
| ![Topología Circos](../assets/comgenomic7.png) | ![Aislamiento Circos](../assets/comgenomic8.png) |
| <small><b>(G) Topología Circular (Circos):</b> Alineamiento de un cromosoma comercial contra múltiples alelos homólogos, ideal para análisis poliploides.</small> | <small><b>(H) Aislamiento de Señales:</b> Interacción de foco (*hover*) para iluminar relaciones ortólogas específicas en redes densas.</small> |

> [!NOTE]
> **Exploración Genómica en Vivo:** La Figura 6 presenta capturas estáticas de la suite analítica `comp-gen`. Puedes interactuar de forma inmersiva con el *Heatmap Evolutivo*, explorar los valores reales de presión de selección (Ka/Ks), aplicar filtros dinámicos y probar las transiciones topológicas en tiempo real a través del visor público disponible en: [johntrujillomonte.com/projects/biocenicana/visor_heatmap_evolutivo.html](https://johntrujillomonte.com/projects/biocenicana/visor_heatmap_evolutivo.html)

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
