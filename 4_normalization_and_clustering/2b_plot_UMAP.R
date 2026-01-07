library(Seurat)
library(ggplot2)
library(stringr)

# cartella dove sono i .rds
rds_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering"
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/umap/"

# cartella di output per le UMAP
dir.create(output_dir, showWarnings = FALSE)

# lista file rds
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Loop su tutti i file
for (file in rds_files) {
  
  # carico l'oggetto
  seu <- readRDS(file)
  
  # estraggo il nome del file senza estensione
  fname <- basename(file)
  sample_name <- str_replace(fname, ".rds", "")
  
  # genero la UMAP
  p <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle(sample_name)
  
  # salvo
  ggsave(
    filename = file.path(output_dir, paste0(sample_name, "_UMAP.png")),
    plot = p,
    width = 7, height = 6, dpi = 300
  )
}


# =====================================================
# Carica oggetto Seurat
# =====================================================
seurat_obj <- readRDS("/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/seurat_res_0.7.rds")
#seurat_obj <- readRDS("G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/3_normalization_and_clustering/seurat_res_0.7.rds")

# Directory di output
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/umap/"
res <- 0.7

# =====================================================
# UMAP generale per cluster
# =====================================================
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, paste0("umap_clusters_res_", res, ".pdf")),
       umap_plot, width = 10, height = 8)

# =====================================================
# UMAP evidenziando singoli campioni
# =====================================================
unique_samples <- unique(seurat_obj@meta.data[["seurat_clusters"]])
cat("Campioni trovati:\n")
print(unique_samples)

for (sample in unique_samples) {
  cells_to_highlight <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[["seurat_clusters"]] == sample]
  
  p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    cells.highlight = cells_to_highlight,
    cols.highlight = "red",
    cols = "lightgray",
    pt.size = 0.8,
    raster = FALSE
  ) +
    ggtitle(paste("Cluster:", sample)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(output_dir, paste0("umap_highlight_", sample, ".pdf")),
    plot = p,
    width = 8,
    height = 6
  )
  cat("✓ Salvato UMAP per campione:", sample, "\n")
}

# UMAP con tutti i campioni colorati
p_all <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "patient_name",
  pt.size = 1.2
) +
  ggtitle("All samples") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "umap_all_samples.pdf"),
  plot = p_all,
  width = 10,
  height = 8
)

# Pannelli separati per campione
p_split <- DimPlot(
  seurat_obj,
  reduction = "umap",
  split.by = "patient_name",
  pt.size = 0.8,
  ncol = 3
) +
  ggtitle("Samples in clustered UMAP") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "umap_split_samples.pdf"),
  plot = p_split,
  width = 15,
  height = 10
)

# =====================================================
# UMAP evidenziando condizioni
# =====================================================
conditions <- unique(seurat_obj@meta.data$condition)
cat("Condition found:\n")
print(conditions)

# Evidenzia una condizione per volta
for (condition in conditions) {
  color <- ifelse(condition == "healed", "#7fbc41", "#fdae61")
  
  cells_to_highlight <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$condition == condition]
  
  p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    cells.highlight = cells_to_highlight,
    cols.highlight = color,
    cols = "lightgray",
    pt.size = 0.3,
    raster = FALSE
  ) +
    ggtitle(paste("Condition:", condition)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(output_dir, paste0("umap_highlight_", condition, ".pdf")),
    plot = p,
    width = 8,
    height = 6
  )
  cat("✓ Salvato UMAP per condizione:", condition, "\n")
}

# Pannelli separati per condizione
p_split_cond <- DimPlot(
  seurat_obj,
  reduction = "umap",
  split.by = "condition",
  pt.size = 0.8,
  label = TRUE,
  ncol = 3
) +
  ggtitle("UMAP by condition") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "umap_condition.pdf"),
  plot = p_split_cond,
  width = 20,
  height = 10
)

# =====================================================
# UMAP con colori personalizzati per condizioni
# =====================================================
custom_colors <- c(
  "healed" = "#7fbc41",
  "not_healed" = "#fdae61"
)

umap_single_plot <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 0.5
)

umap_single_plot_colored <- umap_single_plot +
  scale_color_manual(values = custom_colors) +
  labs(title = "UMAP by condition")

print(umap_single_plot_colored)

ggsave(
  filename = file.path(output_dir, "umap_single_condition_colored.pdf"),
  plot = umap_single_plot_colored,
  width = 10,
  height = 8
)

cat("\n✅ Tutti i campioni e condizioni plottati e salvati!\n")

# =====================================================
# Verifica metadata
# =====================================================
cat("\nMetadata patient_name:\n")
print(table(seurat_obj$patient_name))
