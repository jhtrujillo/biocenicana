# =============================================================================
# Validación del mapa genético con polymapR
# Comparación de grupos de ligamiento y longitudes vs BioJava
#
# Instalar si no está disponible:
#   install.packages("polymapR")
#   install.packages("dplyr")
# =============================================================================

library(polymapR)
library(dplyr)

# ── Configuración ──────────────────────────────────────────────────────────
MAP_PATH <- "genomica_comparativa/mapa_genetico/mapa_biparental_completo.map"
VCF_PATH <- "benchmarks/vcfs/mapa_genetico/AllSamples_variants_geneticmap_pseudochromosomes_standarfilters_minInd_85_single_dosage.vcf"
OUT_DIR  <- "genomica_comparativa/mapa_genetico/validacion"
PLOIDY   <- 10
LOD_MIN  <- 3.0
REC_MAX  <- 0.35

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("=== Validación con polymapR ===\n\n")


# ── 1. Leer el mapa de BioJava ─────────────────────────────────────────────
cat("Cargando mapa de BioJava...\n")
biojava_map <- read.table(MAP_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(biojava_map) <- c("Marker","LG","cM","Chr_Phys","Pos_Phys")
cat(sprintf("  %d marcadores en %d LGs\n\n", nrow(biojava_map), length(unique(biojava_map$LG))))


# ── 2. Leer VCF y construir matriz de dosificaciones ──────────────────────
cat("Leyendo VCF y extrayendo dosificaciones (ACN)...\n")

read_vcf_dosages <- function(vcf_path, target_markers) {
  con <- file(vcf_path, "r")
  samples <- NULL
  dosage_list <- list()

  while (TRUE) {
    line <- readLines(con, n=1, warn=FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##")) next
    if (startsWith(line, "#")) {
      cols <- strsplit(line, "\t")[[1]]
      samples <- cols[10:length(cols)]
      next
    }
    cols <- strsplit(line, "\t")[[1]]
    if (length(cols) < 10) next
    chr <- cols[1]; pos <- cols[2]
    mid <- paste0(chr, "_", pos)
    if (!mid %in% target_markers) next

    fmt_fields <- strsplit(cols[9], ":")[[1]]
    acn_idx <- which(fmt_fields == "ACN")

    genos <- sapply(cols[10:length(cols)], function(g) {
      if (startsWith(g, ".")) return(NA_real_)
      fields <- strsplit(g, ":")[[1]]
      if (length(acn_idx) > 0 && length(fields) >= acn_idx) {
        parts <- strsplit(fields[acn_idx], ",")[[1]]
        val <- suppressWarnings(as.numeric(parts[length(parts)]))
        return(ifelse(is.na(val), NA_real_, val))
      }
      return(NA_real_)
    })
    names(genos) <- samples
    dosage_list[[mid]] <- genos
  }
  close(con)
  cat(sprintf("  %d marcadores leídos del VCF\n", length(dosage_list)))
  return(list(dosages=dosage_list, samples=samples))
}

vcf_data <- read_vcf_dosages(VCF_PATH, biojava_map$Marker)
dosages  <- vcf_data$dosages
samples  <- vcf_data$samples

# Construir matriz: filas = marcadores, columnas = individuos
markers_in_vcf <- intersect(biojava_map$Marker, names(dosages))
cat(sprintf("  %d marcadores con datos de dosificación disponibles\n\n", length(markers_in_vcf)))

dose_matrix <- do.call(rbind, lapply(markers_in_vcf, function(m) dosages[[m]]))
rownames(dose_matrix) <- markers_in_vcf
colnames(dose_matrix) <- samples

# Excluir padres si están presentes (primeras 2 columnas suelen ser padres)
# Los padres aparecen como CC_011940 y CC_01746 según el VCF
parent_cols <- grep("CC_011940|CC_01746", colnames(dose_matrix))
if (length(parent_cols) > 0) {
  offspring_matrix <- dose_matrix[, -parent_cols]
  cat(sprintf("  Excluidos %d padres. %d individuos en la progenie.\n\n",
              length(parent_cols), ncol(offspring_matrix)))
} else {
  offspring_matrix <- dose_matrix
}


# ── 3. Calcular LOD y recombinación con polymapR ───────────────────────────
cat("Calculando LOD y frecuencias de recombinación con polymapR...\n")
cat("(Esto puede tomar varios minutos para muchos marcadores)\n\n")

# polymapR requiere que los datos sean enteros
offspring_int <- round(offspring_matrix)
offspring_int[is.na(offspring_int)] <- -1  # -1 como dato faltante en polymapR

tryCatch({
  # Calcular LOD entre todos los pares de marcadores
  lod_results <- calc_LDmaps(
    dosage_matrix = offspring_int,
    ploidy = PLOIDY,
    ncores = parallel::detectCores() - 1
  )

  cat("Agrupando marcadores en LGs con polymapR...\n")
  lg_assignment <- group_markers(
    LODscores  = lod_results$LOD,
    LG_number  = length(unique(biojava_map$LG)),
    LOD_thres  = LOD_MIN,
    REC_thres  = REC_MAX
  )

  # Guardar asignación de polymapR
  polymapR_lgs <- data.frame(
    Marker = names(lg_assignment),
    LG_polymapR = as.integer(lg_assignment),
    stringsAsFactors = FALSE
  )

  # ── 4. Comparación BioJava vs polymapR ──────────────────────────────────
  cat("\n=== COMPARACIÓN BIOJAVA vs polymapR ===\n\n")

  comparison <- merge(biojava_map[, c("Marker","LG")], polymapR_lgs,
                      by="Marker", all=TRUE)
  comparison$match <- comparison$LG == paste0("LG", comparison$LG_polymapR)

  # Concordancia
  n_compared <- sum(!is.na(comparison$LG) & !is.na(comparison$LG_polymapR))
  n_match    <- sum(comparison$match, na.rm=TRUE)
  cat(sprintf("  Marcadores comparados: %d\n", n_compared))
  cat(sprintf("  Asignación concordante: %d (%.1f%%)\n", n_match, 100*n_match/n_compared))

  # LGs de BioJava vs polymapR
  biojava_n_lgs  <- length(unique(na.omit(comparison$LG)))
  polymapR_n_lgs <- length(unique(na.omit(comparison$LG_polymapR)))
  cat(sprintf("  LGs BioJava:   %d\n", biojava_n_lgs))
  cat(sprintf("  LGs polymapR:  %d\n", polymapR_n_lgs))

  # Guardar comparación
  write.table(comparison, file.path(OUT_DIR, "comparison_biojava_polymapR.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("\n  Tabla de comparación guardada en: %s\n",
              file.path(OUT_DIR, "comparison_biojava_polymapR.tsv")))

}, error = function(e) {
  cat(sprintf("  ⚠ polymapR no pudo completar el análisis: %s\n", e$message))
  cat("  Asegúrate de tener polymapR instalado: install.packages('polymapR')\n")
})


# ── 5. Métricas del mapa de BioJava ────────────────────────────────────────
cat("\n=== MÉTRICAS DEL MAPA BIOJAVA ===\n\n")

lg_stats <- biojava_map %>%
  group_by(LG) %>%
  summarise(
    n_markers    = n(),
    length_cM    = max(cM),
    n_chrs       = n_distinct(Chr_Phys),
    dominant_chr = names(sort(table(Chr_Phys), decreasing=TRUE))[1],
    dom_pct      = round(100 * max(table(Chr_Phys)) / n(), 1),
    span_Mb      = round((max(Pos_Phys) - min(Pos_Phys)) / 1e6, 2),
    cM_per_Mb    = round(ifelse(span_Mb > 0, length_cM / span_Mb, NA), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n_markers))

cat(sprintf("  Total LGs:           %d\n", nrow(lg_stats)))
cat(sprintf("  Total marcadores:    %d\n", sum(lg_stats$n_markers)))
cat(sprintf("  Longitud total:      %.1f cM\n", sum(lg_stats$length_cM)))
cat(sprintf("  Longitud promedio:   %.1f cM/LG\n", mean(lg_stats$length_cM)))
cat(sprintf("  cM/Mb promedio:      %.2f\n", mean(lg_stats$cM_per_Mb, na.rm=TRUE)))
cat(sprintf("  LGs con 1 cromosoma: %d (%.1f%%)\n",
            sum(lg_stats$n_chrs == 1), 100*mean(lg_stats$n_chrs == 1)))

print(lg_stats, n=20)

write.table(lg_stats, file.path(OUT_DIR, "lg_stats_biojava.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── 6. Correlación física-genética por LG ──────────────────────────────────
cat("\n=== CORRELACIÓN SPEARMAN POR LG ===\n\n")

spearman_by_lg <- biojava_map %>%
  group_by(LG) %>%
  filter(n() >= 3) %>%
  summarise(
    n = n(),
    spearman_r = cor(Pos_Phys, cM, method="spearman"),
    dominant_chr = names(sort(table(Chr_Phys), decreasing=TRUE))[1],
    .groups = "drop"
  ) %>%
  mutate(evaluation = case_when(
    spearman_r >= 0.7 ~ "✓ bueno",
    spearman_r >= 0.5 ~ "~ aceptable",
    TRUE              ~ "⚠ revisar"
  )) %>%
  arrange(spearman_r)

cat(sprintf("  LGs con Spearman >= 0.7: %d/%d\n",
            sum(spearman_by_lg$spearman_r >= 0.7), nrow(spearman_by_lg)))
cat(sprintf("  Correlación promedio:    %.3f\n", mean(spearman_by_lg$spearman_r)))

print(spearman_by_lg, n=20)

write.table(spearman_by_lg, file.path(OUT_DIR, "spearman_by_lg_R.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

cat(sprintf("\n✓ Validación con R completada. Resultados en: %s/\n", OUT_DIR))
