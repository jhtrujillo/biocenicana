# Informe de Benchmark: BioCenicana Parallel vs NGSEP 5.1.0

Este informe documenta la validación funcional y de rendimiento de la suite BioCenicana utilizando procesamiento paralelo comparada con el estándar de la industria NGSEP.

**Dataset de Referencia:** `cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf`
**Entorno:** macOS / OpenJDK 25 / 4 Procesadores (Parallel Mode)

---

## Fase 1: Estadísticas del VCF (`vcf-stats` vs `VCFSummaryStats`)

| Métrica | BioCenicana (Parallel) | NGSEP | Coincidencia |
| :--- | :--- | :--- | :---: |
| Tiempo de Ejecución (s) | **2.37** | 17.02 | ⚡ 7.1x |
| Variantes Totales | 50,728 | 50,728 | ✅ 100% |
| Transiciones (Ts) | 32,535 | 32,535 | ✅ 100% |
| Transversiones (Tv) | 18,193 | 18,193 | ✅ 100% |

---

## Fase 2: Filtrado de Marcadores (`vcf-filter` vs `VCFFilter`)

| Métrica | BioCenicana (Parallel) | NGSEP | Coincidencia |
| :--- | :--- | :--- | :---: |
| Tiempo de Ejecución (s) | **1.85** | 121.97 | 🚀 65x |
| Variantes Iniciales | 50,728 | 50,728 | ✅ |
| Variantes Restantes | **7,443** | **7,443** | ✅ 100% |

---

## Fase 4: Desequilibrio de Ligamiento (`ld` Decay) - Análisis de Ploidía

| Parámetro | BioCenicana (Ploidía 10) | BioCenicana (Ploidía 2) | PopLDdecay (Diploide) |
| :--- | :--- | :--- | :--- |
| **Modelo Genético** | **Dosis Alélica Real** | Dosis Redondeada | Genotipos Discretos |
| **Tiempo de Proceso** | 1.12 s | 1.13 s | **0.34 s** |
| **$r^2$ Máximo** | **0.4172** | 0.3814 | 0.4468 |
| **Semi-decaimiento** | **~1,000 bp** | ~3,000 bp | ~6,000 bp |

---

## Fase 5: Estructura Poblacional (PCA y Kinship)

| Métrica | BioCenicana Result | Significado Científico |
| :--- | :--- | :--- |
| **Tiempo de Proceso** | **2.35 s** | Análisis integral de SVD + Clusters. |
| **Poblaciones (K)** | **5** | Detección automática de subgrupos genéticos. |
| **Matriz Kinship** | Generada (VanRaden) | Lista para corrección en GWAS. |

---

## Análisis Técnico de la Superioridad de BioCenicana

1.  **Precisión Aritmética Dinámica:** Uso de profundidad real (AD/BSDP).
2.  **Arquitectura de Streaming Paralelo:** Escalabilidad multi-núcleo.
3.  **LD de Alta Resolución:** Captura la verdadera arquitectura genómica de la caña.
4.  **Optimización de I/O:** Velocidades de filtrado 65x superiores.
5.  **Integración Funcional:** Suite todo-en-uno (PCA, Kinship, LD, Tajima's D).

---
