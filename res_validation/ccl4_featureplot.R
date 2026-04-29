# -----------------------------------------------------------------------------
# Seurat Visualization Pipeline for CCL4
# Purpose: Generate 3 high-visibility FeaturePlots for gene CCL4
# -----------------------------------------------------------------------------

# 1. Setup and Libraries
library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Define paths
seurat_file_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds"
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/UMAPs_CCL4/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data
message("Loading Seurat object... this might take a moment due to file size.")
#seurat_obj <- readRDS(seurat_file_path)

# Verify object structure and gene existence
# We set the default assay to RNA to ensure we access the correct data slots
DefaultAssay(seurat_obj) <- "RNA"

target_gene <- "CCL4"

if (!target_gene %in% rownames(seurat_obj)) {
  stop(paste("Gene", target_gene, "not found in the RNA assay."))
}

# 3. Preparation
# Ensure identity is set to the requested resolution context (though FeaturePlot relies on UMAP coords)
if ("RNA_snn_res.0.7" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- "RNA_snn_res.0.7"
}

# Optimization for visibility:
# For 292k cells, point size must be minimal to avoid overplotting occlusion
pt_size_opt <- 0.1 

message(paste("Generating plots for:", target_gene))

# -----------------------------------------------------------------------------
# OPTION 1: The "Ordered & Clipped" Plot (Red/Grey)
# Strategy: Force expressing cells to top (order=T) and clip top 5% outliers 
# to prevents color scale flattening.
# -----------------------------------------------------------------------------
p1 <- FeaturePlot(
  object = seurat_obj,
  features = target_gene,
  pt.size = pt_size_opt,
  order = TRUE,             # CRITICAL: Plots high expression cells on top of non-expressing
  min.cutoff = "q1",        # Reduces noise at the bottom
  max.cutoff = "q98",       # Clips extreme outliers to brighten the rest
  cols = c("lightgrey", "#D73027"), # High contrast Red
  raster = FALSE            # Keep vector quality
) + 
  ggtitle(paste(target_gene, "Expression")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste0(output_dir, "1_CCL4_Red_Ordered_Clipped.png"), 
       plot = p1, width = 10, height = 8, dpi = 300)

message("Processing complete. Plots saved in: ", output_dir)