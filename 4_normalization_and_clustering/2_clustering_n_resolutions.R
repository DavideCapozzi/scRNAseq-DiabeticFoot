library(Seurat)
library(dplyr)

set.seed(123)

# =========================
# CONFIGURAZIONE
# =========================
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering"
#output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/3_normalization_and_clustering/umap_highlit_condition/"

seurat_obj <- readRDS("/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_pca/seurat_after_normalization_and_pca.rds")
#seurat_obj <- readRDS("G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/3_normalization_and_clustering/seurat_res_0.7.rds")

# Risoluzioni da testare
resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)
#resolutions <- c(0.7)

# =========================
# Step 1: FindNeighbors
# =========================
cat("\n==============================\n")
cat("STEP 1: FIND NEIGHBORS\n")
cat("==============================\n")

# Controlla se il grafo esiste già
if (!"RNA_snn" %in% names(seurat_obj@graphs)) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:6, verbose = TRUE)
  cat("✓ Grafo delle similarità costruito\n")
} else {
  cat("ℹ Grafo già presente, skip FindNeighbors\n")
}

# =========================
# Step 2: Loop su risoluzioni
# =========================
for (res in resolutions) {
  cat("\n==============================\n")
  cat("Clustering con resolution =", res, "\n")
  cat("==============================\n")
  
  # Clustering
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = TRUE)
  n_clusters <- length(unique(seurat_obj$seurat_clusters))
  cat("Clusters identificati:", n_clusters, "\n")
  
  # UMAP
  # Controlla se UMAP esiste già
  if (!"umap" %in% names(seurat_obj@reductions)) {
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:6, verbose = TRUE)
    cat("✓ UMAP calcolato\n")
  } else {
    cat("ℹ UMAP già presente\n")
  }
  
  # --- Salva oggetto Seurat per questa resolution ---
  saveRDS(seurat_obj, file.path(output_dir, paste0("seurat_res_", res, ".rds")))
}

cat("\n✅ Tutte le risoluzioni completate e salvate!\n")
