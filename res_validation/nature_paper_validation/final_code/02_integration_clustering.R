# ==============================================================================
# SCRIPT 02: HARMONY INTEGRATION, CLUSTERING AND CELL TYPE ANNOTATION
# Dataset: GSE165816 - DFU scRNA-seq
#
# WORKFLOW:
#   1. Load pre-processed object (output of Script 01)
#   2. Quantitative selection of optimal PCA dimensions
#   3. Harmony batch correction (by patient / orig.ident)
#   4. FindNeighbors on Harmony embedding  [CORRECTED from original code]
#   5. Multi-resolution clustering
#   6. RunUMAP on Harmony embedding        [CORRECTED from original code]
#   7. Cell type annotation (AddModuleScore, 28 cell types)
#   8. UMAP visualisations
#   9. Save annotated object
#
# PREREQUISITE: Script 01 completed
# ==============================================================================

rm(list = ls())
set.seed(42)

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(scales)
library(ggrepel)

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/nature_paper_validation"
res_dir  <- file.path(base_dir, "res")
plot_dir <- file.path(res_dir, "plots")

dir.create(file.path(plot_dir, "clustering"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plot_dir, "annotation"),  showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# GLOBAL PLOT THEME — black axes, no grid, English labels
# ==============================================================================
theme_dfu <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line         = element_line(color = "black", linewidth = 0.5),
      axis.ticks        = element_line(color = "black"),
      axis.text         = element_text(color = "black"),
      axis.title        = element_text(color = "black"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      strip.background  = element_rect(fill = "grey92", color = "black", linewidth = 0.4),
      strip.text        = element_text(color = "black", face = "bold"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.text       = element_text(color = "black"),
      legend.title      = element_text(color = "black", face = "bold"),
      plot.title        = element_text(face = "bold", color = "black", size = base_size + 1),
      plot.subtitle     = element_text(color = "grey35", size = base_size - 1),
      plot.caption      = element_text(color = "grey50",  size = base_size - 2)
    )
}

# Patch applied on top of Seurat-generated ggplot objects
seurat_theme <- function() {
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black"),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    strip.background  = element_rect(fill = "grey92", color = "black"),
    strip.text        = element_text(color = "black", face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black", face = "bold"),
    plot.title        = element_text(face = "bold", color = "black"),
    plot.subtitle     = element_text(color = "grey35")
  )
}

cond_colors <- c("Healing" = "#2196F3", "Non-Healing" = "#F44336", "Diabetes" = "#4CAF50")

# ==============================================================================
# 1. LOAD PRE-PROCESSED OBJECT
# ==============================================================================
cat("--- Loading pre-processed object ---\n")
seurat_obj <- readRDS(file.path(res_dir, "seurat_obj_pre_integration.rds"))

cat(sprintf("Dimensions: %d genes × %d cells\n", nrow(seurat_obj), ncol(seurat_obj)))
cat("Cells per condition:\n"); print(table(seurat_obj$Condition))
cat("Samples:\n");             print(table(seurat_obj$orig.ident))

# ==============================================================================
# 2. OPTIMAL PCA DIMENSION SELECTION
# ==============================================================================
cat("\n--- Selecting optimal number of PCs ---\n")

pct_var <- seurat_obj[["pca"]]@stdev^2 / sum(seurat_obj[["pca"]]@stdev^2) * 100
cum_var <- cumsum(pct_var)

pc_80pct  <- min(which(cum_var >= 80))              # criterion 1: 80% cumulative variance
pc_plateau <- min(which(diff(pct_var) > -0.1)) + 1  # criterion 2: variance increment plateau
delta2    <- diff(diff(pct_var))
pc_elbow  <- which.max(delta2) + 2                  # criterion 3: second derivative (elbow)

optimal_dims <- as.integer(median(c(pc_80pct, pc_plateau, pc_elbow)))
optimal_dims <- max(15L, min(optimal_dims, 35L))    # bound: [15, 35]

cat(sprintf("  pc_80pct=%d | pc_plateau=%d | pc_elbow=%d → optimal: %d PCs\n",
            pc_80pct, pc_plateau, pc_elbow, optimal_dims))
cat(sprintf("  Cumulative variance at %d PCs: %.1f%%\n", optimal_dims, cum_var[optimal_dims]))

# Annotated elbow plot
elbow_df <- data.frame(PC = seq_along(pct_var), Variance = pct_var, Cumulative = cum_var)

p_elbow <- ggplot(elbow_df, aes(x = PC, y = Variance)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 1.5) +
  geom_vline(xintercept = optimal_dims, linetype = "dashed",
             color = "firebrick", linewidth = 1) +
  annotate("text", x = optimal_dims + 1.2, y = max(pct_var) * 0.75,
           label = paste0("Optimal: PC", optimal_dims),
           color = "firebrick", hjust = 0, size = 3.8, fontface = "bold") +
  labs(title = "Elbow Plot — Optimal PC Selection",
       subtitle = sprintf("%.1f%% cumulative variance at %d PCs",
                          cum_var[optimal_dims], optimal_dims),
       x = "Principal Component", y = "% Variance Explained") +
  theme_dfu(12)

ggsave(file.path(plot_dir, "clustering", "01_ElbowPlot.pdf"), p_elbow, width = 8, height = 5)

# ==============================================================================
# 3. HARMONY BATCH CORRECTION
#    Corrects batch effects across patients (orig.ident) while preserving
#    biological variability across conditions (Healing / Non-Healing / Diabetes)
# ==============================================================================
cat("\n--- Harmony batch correction ---\n")
cat("  Correcting for: orig.ident (patient/sample)\n")
cat("  Preserving: biological condition variance\n")

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  reduction.use     = "pca",
  dims.use      = 1:optimal_dims,
  lambda        = 1,
  nclust        = 100,
  max_iter      = 20,
  verbose       = FALSE
)
cat("  Harmony completed. Embedding 'harmony' available.\n")

# Visualise PCA vs Harmony (mixing by patient)
p_pca_batch  <- DimPlot(seurat_obj, reduction = "pca",
                         group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("PCA — Before Harmony") + NoLegend() + seurat_theme()

p_harm_batch <- DimPlot(seurat_obj, reduction = "harmony",
                         group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("Harmony — After Batch Correction") + seurat_theme()

p_batch <- p_pca_batch | p_harm_batch
ggsave(file.path(plot_dir, "clustering", "02_Harmony_PCA_comparison.pdf"),
       p_batch, width = 14, height = 6)

# ==============================================================================
# 4. FIND NEIGHBORS — ON HARMONY EMBEDDING
#    [CORRECTION: reduction = "harmony" (was commented out in original code)]
# ==============================================================================
cat("\n--- FindNeighbors on Harmony embedding ---\n")

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",   # CORRECTED: uses batch-corrected PCs
  dims      = 1:optimal_dims,
  k.param   = 20,
  verbose   = FALSE
)

# ==============================================================================
# 5. MULTI-RESOLUTION CLUSTERING
# ==============================================================================
cat("\n--- Multi-resolution clustering ---\n")
resolutions <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res,
                              cluster.name = paste0("clusters_res", res),
                              verbose = FALSE)
}
cat(sprintf("  Clustering done at %d resolutions.\n", length(resolutions)))

# ==============================================================================
# 6. UMAP — ON HARMONY EMBEDDING
#    [CORRECTION: reduction = "harmony" (was commented out in original code)]
# ==============================================================================
cat("\n--- RunUMAP on Harmony embedding ---\n")

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction   = "harmony",   # CORRECTED: UMAP reflects batch correction
  dims        = 1:optimal_dims,
  n.neighbors = 30,
  min.dist    = 0.3,
  spread      = 1.0,
  n.epochs    = 500,
  seed.use    = 42,
  verbose     = FALSE
)
cat("  UMAP completed.\n")

# UMAP grid at multiple resolutions (for visual inspection)
p_res_list <- lapply(resolutions, function(res) {
  col_name <- paste0("clusters_res", res)
  n_clust  <- length(unique(seurat_obj@meta.data[[col_name]]))
  DimPlot(seurat_obj, group.by = col_name, pt.size = 0.3,
          label = TRUE, label.size = 3) +
    ggtitle(sprintf("res = %.1f  |  %d clusters", res, n_clust)) +
    NoLegend() +
    seurat_theme() +
    theme(plot.title = element_text(size = 10))
})
p_multi_res <- wrap_plots(p_res_list, ncol = 4)
ggsave(file.path(plot_dir, "clustering", "03_UMAP_multi_resolution.pdf"),
       p_multi_res, width = 20, height = 10)

# ---- Choose resolution ----
# res = 0.4 provides a good balance for DFU skin tissue.
# Adjust 'chosen_resolution' after visual inspection of the plot above.
chosen_resolution <- 0.4
Idents(seurat_obj) <- paste0("clusters_res", chosen_resolution)
seurat_obj$seurat_clusters <- Idents(seurat_obj)
cat(sprintf("  Selected resolution: %.1f → %d clusters\n",
            chosen_resolution, length(levels(Idents(seurat_obj)))))

# ==============================================================================
# 7. CELL TYPE ANNOTATION — AddModuleScore (28 cell types, DFU skin panel)
# ==============================================================================
cat("\n--- Cell type annotation ---\n")

cell_type_markers <- list(
  # EPITHELIUM
  "Keratinocytes_Basal"      = c("KRT5","KRT14","TP63","COL17A1","ITGA6"),
  "Keratinocytes_Suprabasal" = c("KRT1","KRT10","KRTDAP","FLG","LOR","IVL"),
  "Keratinocytes_Cycling"    = c("MKI67","TOP2A","PCNA","KRT5"),
  "Keratinocytes_Wound"      = c("KRT6A","KRT16","KRT17","S100A8"),
  # STROMA
  "Fibroblasts"              = c("DCN","LUM","COL1A1","COL3A1","PDGFRA","FBN1"),
  "Myofibroblasts"           = c("ACTA2","TAGLN","MYL9","POSTN"),
  "Fibroblasts_Reticular"    = c("PRRX1","TWIST2","SFRP2","WNT5A"),
  "Fibroblasts_Papillary"    = c("APCDD1","COL18A1","NTN1"),
  # VASCULAR
  "Endothelial_BV"           = c("PECAM1","VWF","FLT1","CLDN5","ERG"),
  "Endothelial_LV"           = c("LYVE1","PROX1","FLT4","PDPN"),
  "Pericytes"                = c("RGS5","PDGFRB","NOTCH3","MCAM"),
  "Smooth_Muscle"            = c("ACTA2","MYH11","CNN1"),
  # INNATE IMMUNITY
  "Macrophages_M1"           = c("CD68","CD80","IL1B","TNF","CXCL10"),
  "Macrophages_M2"           = c("CD68","CD163","MRC1","TGFB1"),
  "Macrophages_General"      = c("CD68","LYZ","CSF1R","AIF1"),
  "Dendritic_Cells"          = c("FCER1A","CD1C","CLEC10A","HLA-DRA"),
  "Langerhans_Cells"         = c("CD207","CD1A","ITGAE"),
  "Mast_Cells"               = c("TPSAB1","TPSB2","CPA3","KIT"),
  "Neutrophils"              = c("S100A8","S100A9","CXCR2","MMP8"),
  "NK_Cells"                 = c("GNLY","NKG7","GZMB","KLRD1"),
  # ADAPTIVE IMMUNITY
  "T_cells_CD8"              = c("CD3D","CD3E","CD8A","GZMK","PRF1"),
  "T_cells_CD4"              = c("CD3D","CD3E","CD4","IL7R","MAL"),
  "T_cells_Treg"             = c("FOXP3","IL2RA","CTLA4"),
  "B_cells"                  = c("CD79A","MS4A1","CD19"),
  "Plasma_cells"             = c("MZB1","IGHG1","IGKC","SDC1"),
  # OTHER SKIN TYPES
  "Melanocytes"              = c("MLANA","DCT","TYRP1","MITF"),
  "Schwann_Neural"           = c("SOX10","MPZ","S100B","PLP1"),
  "Adipocytes"               = c("ADIPOQ","PLIN1","FABP4")
)

cat(sprintf("  Computing module scores for %d cell types...\n", length(cell_type_markers)))
score_names <- c()

for (ct_name in names(cell_type_markers)) {
  genes_ok  <- intersect(cell_type_markers[[ct_name]], rownames(seurat_obj))
  if (length(genes_ok) >= 2) {
    safe_name <- gsub("[^A-Za-z0-9_]", "_", ct_name)
    seurat_obj <- AddModuleScore(seurat_obj, features = list(genes_ok),
                                  name = paste0(safe_name, "_sc"),
                                  ctrl = 50, seed = 42)
    score_names <- c(score_names, paste0(safe_name, "_sc1"))
  }
}

# Aggregate scores per cluster → assign cell type by max score
score_matrix <- sapply(score_names, function(sc) {
  if (sc %in% colnames(seurat_obj@meta.data))
    tapply(seurat_obj@meta.data[[sc]], seurat_obj$seurat_clusters, mean, na.rm = TRUE)
  else
    rep(NA_real_, length(levels(seurat_obj$seurat_clusters)))
})
colnames(score_matrix) <- gsub("_sc1$", "", colnames(score_matrix))

cell_type_assignment <- apply(score_matrix, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || max(x) < 0.05) return("Unknown")
  names(x)[which.max(x)]
})

# Simplify labels for display
simplify_label <- function(x) gsub("_", " ", gsub("_General|_M1|_M2|_BV|_LV", "", x))
cell_type_display <- simplify_label(cell_type_assignment)
names(cell_type_display) <- names(cell_type_assignment)

Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, cell_type_display)
seurat_obj$CellType          <- as.character(Idents(seurat_obj))
seurat_obj$CellType_detailed <- cell_type_assignment

cat("Annotation completed:\n")
print(sort(table(seurat_obj$CellType), decreasing = TRUE))

# Handle Unknown clusters
unknown_clusters <- names(cell_type_assignment)[cell_type_assignment == "Unknown"]
if (length(unknown_clusters) > 0) {
  cat(sprintf("\n  [INFO] %d Unknown cluster(s): %s\n  Running FindMarkers for manual inspection...\n",
              length(unknown_clusters), paste(unknown_clusters, collapse = ", ")))
  Idents(seurat_obj) <- "seurat_clusters"
  tryCatch({
    markers_unknown <- FindMarkers(seurat_obj, ident.1 = unknown_clusters,
                                   min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE)
    write.csv(markers_unknown, file.path(res_dir, "markers_unknown_clusters.csv"))
    cat("  Markers saved: markers_unknown_clusters.csv\n")
  }, error = function(e) cat("  FindMarkers for Unknown failed:", e$message, "\n"))
}
Idents(seurat_obj) <- "CellType"

# ==============================================================================
# 8. UMAP VISUALISATIONS
# ==============================================================================
cat("\n--- Generating UMAP plots ---\n")

n_ct <- length(unique(seurat_obj$CellType))
ct_colors <- scales::hue_pal()(n_ct)
names(ct_colors) <- sort(unique(seurat_obj$CellType))

# UMAP — cell type annotation
p_umap_ct <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "CellType",
  label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.3
) +
  scale_color_manual(values = ct_colors) +
  labs(title = "Cell Type Annotation — DFU Skin (GSE165816)",
       subtitle = sprintf("n = %s cells | %d types | Harmony-integrated UMAP",
                          format(ncol(seurat_obj), big.mark = ","), n_ct),
       x = "UMAP 1", y = "UMAP 2") +
  seurat_theme() + NoLegend()

# UMAP — condition
p_umap_cond <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "Condition",
  pt.size = 0.3, cols = cond_colors
) +
  labs(title = "UMAP by Condition", x = "UMAP 1", y = "UMAP 2") +
  seurat_theme()

# UMAP — patient (check Harmony mixing)
p_umap_patient <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = 0.3
) +
  labs(title = "UMAP by Patient (post-Harmony)", x = "UMAP 1", y = "UMAP 2") +
  seurat_theme()

# UMAP — split by condition (3 panels side by side)
p_umap_split <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "CellType",
  split.by = "Condition", pt.size = 0.3, label = FALSE
) +
  scale_color_manual(values = ct_colors) +
  labs(title = "Cell Types split by Condition", x = "UMAP 1", y = "UMAP 2") +
  seurat_theme() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

# Composite panel
p_panel <- (p_umap_ct | (p_umap_cond / p_umap_patient)) +
  plot_annotation(title = "scRNA-seq DFU — GSE165816 | Harmony-integrated UMAP")
ggsave(file.path(plot_dir, "annotation", "01_UMAP_overview.pdf"),    p_panel,        width = 18, height = 8)
ggsave(file.path(plot_dir, "annotation", "02_UMAP_CellType.pdf"),    p_umap_ct,      width = 10, height = 8)
ggsave(file.path(plot_dir, "annotation", "03_UMAP_split_condition.pdf"), p_umap_split, width = 18, height = 7)

# Dot plot: canonical markers (annotation validation)
canonical_markers <- c("KRT5","KRT14","KRT1","KRT10","DCN","COL1A1",
                        "PECAM1","VWF","LYVE1","RGS5","CD68","LYZ",
                        "TPSAB1","CD3D","CD8A","CD4","GNLY","NKG7",
                        "CD79A","MZB1","MLANA","SOX10","MKI67")
genes_canonical <- canonical_markers[canonical_markers %in% rownames(seurat_obj)]

p_dotplot <- DotPlot(
  seurat_obj, features = genes_canonical, group.by = "CellType",
  dot.scale = 6, col.min = -1, col.max = 2
) +
  coord_flip() +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0,
                         name = "Avg. expression\n(scaled)") +
  labs(title = "Canonical Markers — Annotation Validation",
       x = "Gene", y = "Cell Type") +
  seurat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "annotation", "04_DotPlot_canonical_markers.pdf"),
       p_dotplot, width = 16, height = 10)

# Cell type proportions per condition
cell_prop <- seurat_obj@meta.data %>%
  group_by(Condition, CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))

p_prop <- ggplot(cell_prop, aes(x = Condition, y = prop, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = ct_colors) +
  labs(title = "Cell Type Proportions per Condition",
       y = "Proportion", x = "Condition") +
  theme_dfu(12)
ggsave(file.path(plot_dir, "annotation", "05_CellType_proportions.pdf"),
       p_prop, width = 10, height = 7)

# ==============================================================================
# 9. SAVE ANNOTATED OBJECT
# ==============================================================================
rds_path <- file.path(res_dir, "seurat_obj_annotated.rds")
saveRDS(seurat_obj, rds_path)

cat("\n=== INTEGRATION SUMMARY ===\n")
cat(sprintf("  Optimal PCA dims (Harmony)  : %d\n", optimal_dims))
cat(sprintf("  Clustering resolution       : %.1f\n", chosen_resolution))
cat(sprintf("  Total clusters              : %d\n", length(unique(seurat_obj$seurat_clusters))))
cat(sprintf("  Cell types annotated        : %d\n", n_ct))
cat(sprintf("  Total cells                 : %d\n", ncol(seurat_obj)))
cat(sprintf("  Object saved to             : %s\n", rds_path))
cat("\nProceed with: 03_ccl4_analysis.R\n")
