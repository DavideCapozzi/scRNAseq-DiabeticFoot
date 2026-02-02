# ========================================
# PREPROCESSING SEURAT - ORDINE CORRETTO
# ========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)

set.seed(123)

# ========================================
# CONFIGURAZIONE DIRECTORY
# ========================================
input_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/3_filtering/"
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/3_filtering/"

# input_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/3_filtering/"
# output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_pca/"

# Crea directory output se non esiste
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ========================================
# CARICAMENTO DATI FILTRATI
# ========================================
cat("\n=== LOADING FILTERED DATA ===\n")
input_file <- file.path(input_dir, "merged_seurat_after_QC.rds")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

seurat_obj <- readRDS(input_file)
cat("Loaded filtered Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")

# ========================================
# STEP 1: NORMALIZZAZIONE
# ========================================
cat("\n" , paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("STEP 1: NORMALIZZAZIONE\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("Eseguendo NormalizeData...\n")
seurat_obj <- NormalizeData(
  object = seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)
cat("✓ Normalizzazione completata\n")

# ========================================
# STEP 2: SELEZIONE FEATURES VARIABILI
# ========================================
cat("\n" , paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("STEP 2: SELEZIONE FEATURES VARIABILI\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("Eseguendo FindVariableFeatures...\n")
seurat_obj <- FindVariableFeatures(
  object = seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

# Identifica le top variable features
top10 <- head(VariableFeatures(seurat_obj), 10)
cat("✓ Variable features identificate:", length(VariableFeatures(seurat_obj)), "\n")
cat("Top 10 variable features:\n")
print(top10)

# Plot variable features
var_plot <- VariableFeaturePlot(seurat_obj)
var_plot <- LabelPoints(plot = var_plot, points = top10, repel = TRUE)
ggsave(file.path(output_dir, "variable_features_plot.pdf"), var_plot, width = 10, height = 6)
cat("✓ Variable features plot salvato\n")

# ========================================
# STEP 3: SCALING DEI DATI
# ========================================
cat("\n" , paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("STEP 3: SCALING DEI DATI\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("Eseguendo ScaleData...\n")
seurat_obj <- ScaleData(
  object = seurat_obj,
  features = VariableFeatures(seurat_obj),  # Scala solo le variable features
  verbose = TRUE
)
cat("✓ Scaling completato\n")

# ========================================
# STEP 4: ANALISI PCA
# ========================================
cat("\n" , paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("STEP 4: ANALISI PCA\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("Eseguendo RunPCA...\n")
seurat_obj <- RunPCA(
  object = seurat_obj,
  features = VariableFeatures(seurat_obj),  # PCA solo sulle variable features
  npcs = 50,
  verbose = TRUE
)
cat("✓ PCA completata con 50 componenti principali\n")

# Numero di PC da utilizzare per downstream analysis
n_pcs <- 6

# Calcola e stampa la varianza spiegata
pca_variance <- seurat_obj@reductions$pca@stdev^2
total_variance <- sum(pca_variance)
variance_explained <- sum(pca_variance[1:n_pcs]) / total_variance * 100

cat("\n📊 VARIANZA SPIEGATA:\n")
cat(sprintf("  - Prime %d PC: %.2f%%\n", n_pcs, variance_explained))
cat(sprintf("  - Varianza cumulativa per PC:\n"))
for (i in 1:n_pcs) {
  cum_var <- sum(pca_variance[1:i]) / total_variance * 100
  cat(sprintf("    PC 1-%d: %.2f%%\n", i, cum_var))
}

# ========================

# Visualizzazione PCA
cat("Generando visualizzazioni PCA...\n")

# Elbow plot per determinare dimensionalità
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(output_dir, "pca_elbow_plot.pdf"), elbow_plot, width = 8, height = 6)

# PCA plot colored by sample (se disponibile)
if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
  pca_sample_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident")
  ggsave(file.path(output_dir, "pca_by_sample.pdf"), pca_sample_plot, width = 10, height = 8)
}

# Heatmap dei primi PC
pca_heatmap <- DimHeatmap(seurat_obj, dims = 1:9, cells = 500, balanced = TRUE)
ggsave(file.path(output_dir, "pca_heatmap.pdf"), pca_heatmap, width = 12, height = 10)

cat("✓ Plot PCA salvati\n")

output_file <- file.path(output_dir, "seurat_after_normalization_and_pca.rds")
saveRDS(seurat_obj, output_file)
cat("✓ Oggetto Seurat salvato:", output_file, "\n")
