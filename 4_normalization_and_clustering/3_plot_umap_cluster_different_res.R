library(Seurat)
library(ggplot2)
library(stringr)

# Directory dove sono i .rds
#rds_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering"
rds_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"


# Directory di output
output_dir <- file.path(rds_dir, "umap_all_resolutions")
dir.create(output_dir, showWarnings = FALSE)

# Lista file .rds
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

for (file in rds_files) {
  
  seu <- readRDS(file)
  fname <- basename(file)
  sample_name <- str_replace(fname, ".rds", "")
  
  message("\n=== Elaboro risoluzione: ", sample_name, " ===")
  
  # Sottocartella per questa risoluzione
  out_subdir <- file.path(output_dir, sample_name)
  dir.create(out_subdir, showWarnings = FALSE)
  
  # =====================================================
  # 1) UMAP NORMALE
  # =====================================================
  p_norm <- DimPlot(seu, reduction = "umap", pt.size = 0.6) +
    ggtitle(paste("UMAP -", sample_name))
  
  ggsave(file.path(out_subdir, paste0(sample_name, "_UMAP_normal.png")),
         p_norm, width = 7, height = 6, dpi = 300)
  
  # =====================================================
  # 2) UMAP NORMALE CON LABEL
  # =====================================================
  p_label <- DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.6) +
    ggtitle(paste("UMAP with labels -", sample_name))
  
  ggsave(file.path(out_subdir, paste0(sample_name, "_UMAP_labeled.png")),
         p_label, width = 7, height = 6, dpi = 300)
  
  # =====================================================
  # 3) UMAP CLUSTER PER CLUSTER (rosso vs grigio)
  # =====================================================
  clusters <- sort(unique(seu$seurat_clusters))
  message("Cluster trovati: ", paste(clusters, collapse = ", "))
  
  for (cl in clusters) {
    
    cells_highlight <- rownames(seu@meta.data)[seu$seurat_clusters == cl]
    
    p_cl <- DimPlot(
      seu,
      reduction = "umap",
      cells.highlight = cells_highlight,
      cols.highlight = "red",
      cols = "lightgray",
      pt.size = 0.6
    ) +
      ggtitle(paste(sample_name, "- Cluster", cl)) +
      theme_minimal()
    
    ggsave(
      filename = file.path(out_subdir, paste0(sample_name, "_highlight_cluster_", cl, ".png")),
      plot = p_cl,
      width = 7,
      height = 6, dpi = 300
    )
  }
  
  # =====================================================
  # 4) UMAP PER CAMPIONE (COLORE = patient_name)
  # =====================================================
  if ("patient_name" %in% colnames(seu@meta.data)) {
    
    p_samples <- DimPlot(
      seu,
      reduction = "umap",
      group.by = "patient_name",
      pt.size = 0.7
    ) +
      ggtitle(paste("Samples colored -", sample_name)) +
      theme_minimal()
    
    ggsave(
      filename = file.path(out_subdir, paste0(sample_name, "_UMAP_by_sample.png")),
      plot = p_samples,
      width = 8,
      height = 7,
      dpi = 300
    )
    
    message(" ✓ Salvata UMAP colorata per campione")
  } else {
    message(" ⚠ patient_name non trovato nei metadata.")
  }
  
  # =====================================================
  # 5) UMAP SPLIT-BY-SAMPLE (pannelli)
  # =====================================================
  if ("patient_name" %in% colnames(seu@meta.data)) {
    
    p_split_samples <- DimPlot(
      seu,
      reduction = "umap",
      split.by = "patient_name",
      pt.size = 0.6,
      ncol = 3
    ) +
      ggtitle(paste("Split by sample -", sample_name)) +
      theme_minimal()
    
    ggsave(
      filename = file.path(out_subdir, paste0(sample_name, "_UMAP_split_samples.png")),
      plot = p_split_samples,
      width = 15,
      height = 10,
      dpi = 300
    )
    
    message(" ✓ Salvata UMAP split-by-sample")
  }
}

message("\n🎉 TUTTE le UMAP sono state generate, incluse quelle per campione!")

