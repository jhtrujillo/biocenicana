# Biocenicana

Herramientas bioinformáticas para genómica de poblaciones y estudios de huella molecular en organismos poliploides (como caña de azúcar). Construido con Java, Maven y [PicoCLI](https://picocli.info/).

## Características Principales

*   **`allele-dosage`**: Genera rápidamente matrices de dosis alélicas a partir de archivos VCF masivos, sin problemas de sobrecarga de memoria (`OutOfMemoryError`). Soporta detección automática de VCFs provenientes de NGSEP, GATK y FreeBayes, permitiendo aplicar modelos de imputación dinámicos al vuelo.
*   **Lectura VCF Optimizada**: Implementa un lector línea por línea en `VcfFastReader.java` que garantiza un uso mínimo de RAM, ideal para procesar genomas altamente repetitivos/poliploides.

## Requisitos

*   **Java 11** o superior.
*   **Apache Maven 3.6+** (para compilación).

*(Nota: La dependencia de la librería NGSEPcore viene incluida localmente en la carpeta `/lib` del proyecto).*

## Compilación

Para compilar el proyecto y generar el archivo JAR ejecutable con todas las dependencias, ejecuta el siguiente comando en la raíz del proyecto:

```bash
mvn clean package -DskipTests
```

El proceso generará en la carpeta `target/` un jar empaquetado y listo para correr: `biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar`.

## Uso

El archivo empaquetado es un ejecutable que engloba un toolkit de comandos (`biocenicana`). Te recomendamos crear un alias para invocarlo más fácilmente, o correrlo directamente usando `java -jar`.

### Ver la ayuda general
```bash
java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar --help
```

### Cálculo de Matriz de Dosis Alélica (`allele-dosage`)

Calcula la matriz de genotipos convertidos a dosis alélicas para una ploidía dada, extraída desde la columna FORMAT del VCF.

```bash
java -jar target/biocenicana-1.0-SNAPSHOT-jar-with-dependencies.jar allele-dosage \
  --vcf ruta/al/archivo.vcf \
  --ploidy 10 \
  --impute bsdp-mode
```

#### Opciones útiles:
*   `-p, --ploidy`: Define la ploidía del organismo. (Ej: `10` para caña). Por defecto es `2`.
*   `-i, --impute`: Método de imputación de valores faltantes (missing).
    *   `bsdp`: Deja el valor faltante en `-1`.
    *   `mode`: Imputa según la moda de dicho SNP.
    *   `mean`: Imputa según la media de dicho SNP.
    *   `bsdp-mode` / `bsdp-mean`: Busca conteos primero, y los `-1` restantes los imputa en moda o media.
*   `-c, --caller`: Define de qué software de Variant Calling proviene el VCF. Opciones: `auto`, `ngsep`, `gatk`, `freebayes`. Por defecto es `auto`.

## Estructura del Código

Este repositorio fue recientemente modernizado para abandonar scripts monolíticos en favor de una CLI limpia basada en métodos de clase:
*   `src/main/java/huellamolecular/cli/`: Contiene los comandos individuales de consola invocados por PicoCLI.
*   `src/main/java/huellamolecular/io/`: Parsers eficientes como `VcfFastReader` para evitar caídas de RAM.
*   `src/main/java/huellamolecular/`: Clases lógicas analíticas legacy actualizadas como `GeneDosis.java`.
*   `lib/`: JARs dependencias críticas (`NGSEPcore`).

---
Proyecto gestionado y mantenido por [jhtrujillo](https://github.com/jhtrujillo) y el equipo de investigación de Cenicaña.
