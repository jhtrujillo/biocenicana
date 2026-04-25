# Informe de Validación Paso a Paso: BioCenicana Parallel vs NGSEP

Este informe documenta la validación detallada de la suite BioCenicana utilizando procesamiento paralelo (4 hilos) comparada con NGSEP 5.1.0 sobre el dataset original de producción.

**Dataset de Referencia:** `cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf`
**Tamaño:** 50,728 variantes / 220 muestras
**Entorno:** macOS / OpenJDK 25 / 4 Cores asignados

---

## Fase 1: Estadísticas del VCF (`vcf-stats` vs `VCFSummaryStats`)

Validación de la exactitud del motor de lectura y diagnóstico genético.

### Comandos Ejecutados

**BioCenicana (4t):**
```bash
time java -jar target/biocenicana-1.0.jar vcf-stats \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -o stats_biocenicana_50k -p 10 -t 4
```

**NGSEP:**
```bash
time java -jar benchmarks/NGSEPcore_5.1.0.jar VCFSummaryStats \
  -i benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -o stats_ngsep_50k.txt
```

### Resultados de Rendimiento y Consistencia (Fase 1)

| Métrica | BioCenicana (4t) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :---: |
| Tiempo de Ejecución (s) | **2.37** | 17.02 | ⚡ 7.1x |
| Variantes Totales | 50,728 | 50,728 | ✅ 100% |
| Transiciones (Ts) | 32,535 | 32,535 | ✅ 100% |
| Transversiones (Tv) | 18,193 | 18,193 | ✅ 100% |
| Ts/Tv Ratio | 1.788 | 1.79 | ✅ 100% |

---

## Fase 2: Filtrado de Marcadores (`vcf-filter` vs `VCFFilter`)

Se aplicó el filtro de calidad: **Mínimo 200 muestras genotipadas** y **Bialélicos**.

### Comandos Ejecutados

**BioCenicana (4t + Strict):**
```bash
time java -jar target/biocenicana-1.0.jar vcf-filter \
  -v benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -o filtered_50k_strict_biocenicana.vcf \
  -m 200 -b -p 10 -t 4 -s
```

**NGSEP:**
```bash
time java -jar benchmarks/NGSEPcore_5.1.0.jar VCFFilter \
  -i benchmarks/sugarcane/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf \
  -o filtered_50k_ngsep.vcf \
  -m 200 -s
```

### Resultados de Rendimiento y Consistencia (Fase 2)

| Métrica | BioCenicana (4t + Strict) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :---: |
| Tiempo de Ejecución (s) | **1.85** | 121.97 | 🚀 65x |
| Variantes Iniciales | 50,728 | 50,728 | ✅ |
| Variantes Restantes | **7,443** | **7,443** | ✅ 100% |

---

## Fase 3: Distancia Genética (`genetic-distance` vs `VCFDistanceMatrixCalculator`)

Cálculo de la matriz de distancias genéticas para las 220 muestras utilizando 7,443 SNPs.

### Comandos Ejecutados

**BioCenicana (4t):**
```bash
time java -jar target/biocenicana-1.0.jar genetic-distance \
  -v filtered_50k_strict_biocenicana.vcf \
  -o dist_biocenicana_50k.tsv \
  -p 10 -t 4
```

**NGSEP:**
```bash
time java -jar benchmarks/NGSEPcore_5.1.0.jar VCFDistanceMatrixCalculator \
  -i filtered_50k_ngsep.vcf \
  -o dist_ngsep_50k.tsv \
  -s 3 -p 10
```

### Resultados de Rendimiento y Consistencia (Fase 3)

| Métrica | BioCenicana (4t) | NGSEP 5.1.0 | Coincidencia |
| :--- | :--- | :--- | :---: |
| Tiempo de Ejecución (s) | **1.21** | 4.72 | ⚡ 3.9x |
| Tamaño de Matriz | 220 x 220 | 220 x 220 | ✅ |
| Consistencia de Ranking | **Idéntica** | **Idéntica** | ✅ 100% |

---
