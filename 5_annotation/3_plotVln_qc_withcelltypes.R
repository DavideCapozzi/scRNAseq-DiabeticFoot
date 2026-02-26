# ============================================================
# Violin plots - QC metrics by cell type
# Requires: seurat_obj already loaded in environment
# ============================================================

library(Seurat)
library(patchwork)
library(ggplot2)

# ---- USER SETTINGS -----------------------------------------
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/2_pre_filtering/qc_violin_plots_withcelltypes"   
# ------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Colori per i 19 cluster
n_clusters    <- length(levels(seurat_obj@active.ident))
cluster_cols  <- setNames(
  scales::hue_pal()(n_clusters),
  levels(seurat_obj@active.ident)
)

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo")

plots <- lapply(features, function(feat) {
  VlnPlot(
    object   = seurat_obj,
    features = feat,
    pt.size  = 0,
    cols     = cluster_cols
  ) +
    labs(title = feat, x = NULL, y = feat) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y    = element_blank(),
      legend.position = "none"
    )
})

final_plot <- wrap_plots(plots, ncol = 1)

# Salva
out_file <- file.path(output_dir, "violin_qc_by_celltype.pdf")
ggsave(
  filename = out_file,
  plot     = final_plot,
  width    = 16,
  height   = 14,
  device   = "pdf"
)

message("Plot salvato in: ", out_file)