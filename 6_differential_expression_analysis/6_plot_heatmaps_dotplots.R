# ==============================================================================
# SCRIPT: 6_plot_targetgenes_heatmaps.R
# PURPOSE: Generate ordered, publication-ready visualizations comparing 
#          'healed' vs 'not_healed'. Includes Intra-Cluster Z-Score Heatmap 
#          and a Raw Mean Expression Heatmap variant.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. DEPENDENCIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 1. CONFIGURATION BLOCK
# ------------------------------------------------------------------------------
CONFIG <- list(
  # File paths
  input_seurat   = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds",
  output_dir     = "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/violinplots_our_sc_data/global_comparisons/new_differential_plots/",
  
  # Target genes
  target_genes   = c("CCL4", "CCR5"),
  
  # Heatmap aesthetics
  # Divergent scale for Z-Score (Centered at 0)
  zscore_colors  = c("low" = "#2166AC", "mid" = "#F7F7F7", "high" = "#B2182B"),
  # Sequential scale for Raw Expression (Starts at 0)
  raw_colors     = c("low" = "#F7F7F7", "high" = "#B2182B"),
  
  # Biological ordering of cell lineages (T cells -> NK -> B -> Myeloid)
  cell_lineage_order = c(
    "CD4 T cells (naive/CM)",
    "CD4 T cells (CD4 memory/activated Th2 cells)",
    "CD4 T cells (activated CD4 memory T cells)",
    "CD8 T cells (EM)",
    "CD4/CD8 T cells (memory T cells)",
    "T cells (cytotoxic gamma delta)",
    "T cells (NKT-like gamma delta)",
    "NK cells (CD16)",
    "B cells (naive/transitional)",
    "B cells (activated/pre-plasmablasts)",
    "Classical monocytes",
    "Intermediate monocytes",
    "Nonclassical monocytes",
    "pDCs & cDC2"
  ),
  
  # Plot dimensions
  plot_width     = 16, 
  plot_height    = 6
)

# ------------------------------------------------------------------------------
# 2. WRAPPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Generate an Intra-Cluster Z-Score Expression Heatmap
#' 
#' @param hm_data Pre-fetched long dataframe.
#' @param genes Target genes.
#' @param cell_order Ordered factors for cell types.
#' @return ggplot object.
generate_zscore_heatmap <- function(hm_data, genes, cell_order) {
  
  sc_scaled <- hm_data %>%
    group_by(ident, Gene) %>%
    mutate(
      Scaled_Expr = if (sd(Expression, na.rm = TRUE) > 0) {
        as.vector(scale(Expression))
      } else {
        0 
      }
    ) %>%
    ungroup()
  
  avg_std <- sc_scaled %>%
    group_by(ident, condition, Gene) %>%
    summarise(Mean_Z_Score = mean(Scaled_Expr, na.rm = TRUE), .groups = 'drop')
  
  avg_std$Gene <- factor(avg_std$Gene, levels = rev(genes))
  avg_std$ident <- factor(avg_std$ident, levels = cell_order)
  avg_std$condition <- factor(avg_std$condition, levels = c("healed", "not_healed"))
  
  avg_std <- avg_std %>% filter(!is.na(ident))
  
  max_z <- max(abs(avg_std$Mean_Z_Score), na.rm = TRUE)
  limit_range <- c(-max_z, max_z)
  
  p_heat <- ggplot(avg_std, aes(x = condition, y = Gene, fill = Mean_Z_Score)) +
    geom_tile(color = "white", size = 0.5) +
    facet_grid(~ ident, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = CONFIG$zscore_colors["low"], 
      mid = CONFIG$zscore_colors["mid"], 
      high = CONFIG$zscore_colors["high"], 
      midpoint = 0, 
      limits = limit_range,
      name = "Intra-Cluster\nZ-Score"
    ) +
    # UI FIX: expand = c(0,0) removes native ggplot padding around tiles
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      # UI FIX: Strip all spacing between top text and the tiles
      axis.title.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.ticks.length.x.top = unit(0, "pt"),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, face = "plain", size = 10, margin = margin(b = 2)),
      
      axis.text.y = element_text(face = "bold.italic", size = 12),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.placement = "outside",
      strip.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 10),
      strip.background = element_blank(),
      legend.position = "right"
    )
  
  return(p_heat)
}

#' Generate Raw Mean Expression Heatmap
#' Plots the true log-normalized expression values without scaling.
#' Uses a sequential color scale.
#' 
#' @param hm_data Pre-fetched long dataframe.
#' @param genes Target genes.
#' @param cell_order Ordered factors for cell types.
#' @return ggplot object.
generate_raw_expression_heatmap <- function(hm_data, genes, cell_order) {
  
  # Calculate mean expression per group directly
  avg_expr <- hm_data %>%
    group_by(ident, condition, Gene) %>%
    summarise(Mean_Expr = mean(Expression, na.rm = TRUE), .groups = 'drop')
  
  avg_expr$Gene <- factor(avg_expr$Gene, levels = rev(genes))
  avg_expr$ident <- factor(avg_expr$ident, levels = cell_order)
  avg_expr$condition <- factor(avg_expr$condition, levels = c("healed", "not_healed"))
  
  avg_expr <- avg_expr %>% filter(!is.na(ident))
  
  p_raw <- ggplot(avg_expr, aes(x = condition, y = Gene, fill = Mean_Expr)) +
    geom_tile(color = "white", size = 0.5) +
    facet_grid(~ ident, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradient(
      low = CONFIG$raw_colors["low"], 
      high = CONFIG$raw_colors["high"], 
      name = "Mean Expression\n(Log-norm)"
    ) +
    # UI FIX: expand = c(0,0) removes native ggplot padding around tiles
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      # UI FIX: Strip all spacing between top text and the tiles
      axis.title.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.ticks.length.x.top = unit(0, "pt"),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, face = "plain", size = 10, margin = margin(b = 2)),
      
      axis.text.y = element_text(face = "bold.italic", size = 12),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.placement = "outside",
      strip.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 10),
      strip.background = element_blank(),
      legend.position = "right"
    )
  
  return(p_raw)
}

# ------------------------------------------------------------------------------
# 3. MAIN EXECUTION BLOCK
# ------------------------------------------------------------------------------
tryCatch({
  
  message("Initializing workspace and validating paths...")
  if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE)
  }
  
  message("Loading Seurat object...")
  #seurat_obj <- readRDS(CONFIG$input_seurat) # Uncomment in your real run
  
  message("Sanitizing metadata...")
  seurat_obj$condition <- tolower(seurat_obj$condition)
  
  current_idents <- as.character(Idents(seurat_obj))
  clean_idents <- str_remove_all(current_idents, "^\\d+\\.\\s*")
  Idents(seurat_obj) <- factor(clean_idents, levels = CONFIG$cell_lineage_order)
  
  valid_genes <- intersect(CONFIG$target_genes, rownames(seurat_obj))
  
  message("Fetching data matrix...")
  hm_data_raw <- FetchData(seurat_obj, vars = c(valid_genes, "ident", "condition"), slot = "data")
  hm_data_long <- hm_data_raw %>%
    pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression")
  
  message("Generating Intra-Cluster Z-Score Heatmap...")
  plot_zscore <- generate_zscore_heatmap(hm_data_long, valid_genes, CONFIG$cell_lineage_order)
  
  message("Generating Raw Mean Expression Heatmap...")
  plot_raw <- generate_raw_expression_heatmap(hm_data_long, valid_genes, CONFIG$cell_lineage_order)
  
  message("Saving plots...")
  ggsave(
    filename = file.path(CONFIG$output_dir, "Target_Genes_IntraCluster_ZScore_Ordered.pdf"),
    plot = plot_zscore,
    width = CONFIG$plot_width,
    height = CONFIG$plot_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(CONFIG$output_dir, "Target_Genes_RawExpression_Ordered.pdf"),
    plot = plot_raw,
    width = CONFIG$plot_width,
    height = CONFIG$plot_height,
    dpi = 300
  )
  
  message("Execution completed successfully. Plots are saved in: ", CONFIG$output_dir)
  
}, error = function(e) {
  message("\n[CRITICAL ERROR] Script execution halted.")
  message("Details: ", e$message)
})