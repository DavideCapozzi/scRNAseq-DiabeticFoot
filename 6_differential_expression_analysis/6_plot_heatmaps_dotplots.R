# ==============================================================================
# SCRIPT: 6_plot_heatmaps_dotplots.R
# PURPOSE: Publication-ready visualisations comparing 'healed' vs 'not_healed'.
#
#   For each gene pair (CCL4/CCR5, CCR1/CCR8) the following PDFs are produced:
#     1.  <TAG>_ZScore_Heatmap_Horizontal.pdf
#         Horizontal heatmap — rows = genes, columns = conditions,
#         facets = cell types (original layout).
#
#     2.  <TAG>_ZScore_Heatmap_Vertical.pdf
#         Vertical heatmap — one narrow panel per gene, side-by-side.
#         rows = cell types, columns = healed / not_healed.
#         Tiles annotated with significance asterisks (Wilcoxon rank-sum test).
#
#     3.  <TAG>_DotPlot.pdf
#         Custom dot plot — one panel per gene.
#         rows = cell types, x = condition (healed | not_healed).
#         Dot SIZE = % cells expressing the gene (scaled per gene [1,8]).
#         Dot COLOUR = mean normalised expression.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. DEPENDENCIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)   # side-by-side panel assembly
})


# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------
CONFIG <- list(
  
  # ---- Paths ----------------------------------------------------------------
  input_seurat = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds",
  output_dir   = "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/heatmaps_dotplots/",
  
  # ---- Gene pairs: one full set of PDFs per pair ---------------------------
  gene_pairs = list(
    c("CCL4", "CCR5"),
    c("CCR1", "CCR8")
  ),
  
  # ---- Colour scales --------------------------------------------------------
  zscore_colors = c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B"),
  expr_colors   = c(low = "#F0F0F0", high = "#B2182B"),   # DotPlot expression
  
  # ---- Significance thresholds for Wilcoxon test (used in vertical heatmap)
  sig_thresholds = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
  
  # ---- Biological cell-type order (T -> NK -> B -> Myeloid) ----------------
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
  
  # ---- Output dimensions (inches) ------------------------------------------
  heatmap_h_width  = 18,
  heatmap_h_height = 4,
  heatmap_v_width  = 9,    # narrow: two slim panels side-by-side
  heatmap_v_height = 11,
  dotplot_width    = 14,
  dotplot_height   = 11
)


# ==============================================================================
# 2. UTILITY FUNCTIONS
# ==============================================================================

# Build a clean file-name prefix from a gene vector (e.g. "CCL4_CCR5")
genes_to_tag <- function(genes) paste(toupper(genes), collapse = "_")

# Add intra-cluster Z-score column to a long expression data frame.
# Groups by (ident, Gene); adds 'Scaled_Expr'.
add_zscore <- function(df) {
  df %>%
    group_by(ident, Gene) %>%
    mutate(
      Scaled_Expr = if (sd(Expression, na.rm = TRUE) > 0) {
        as.vector(scale(Expression))
      } else {
        0
      }
    ) %>%
    ungroup()
}

# Convert a numeric p-value to a significance asterisk string.
# Returns "***", "**", "*", or "" (empty string = not significant).
pval_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# Run Wilcoxon rank-sum tests (healed vs not_healed) for every (cell type, gene)
# combination present in a long expression data frame.
# Returns a data frame with columns: ident, Gene, p_value, stars.
compute_wilcoxon <- function(long_df) {
  long_df %>%
    group_by(ident, Gene) %>%
    summarise(
      p_value = {
        grp_h  <- Expression[condition == "healed"]
        grp_nh <- Expression[condition == "not_healed"]
        # Require at least 3 cells per group; otherwise return NA
        if (length(grp_h) >= 3 && length(grp_nh) >= 3) {
          tryCatch(
            wilcox.test(grp_h, grp_nh, exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(stars = pval_to_stars(p_value))
}


# ==============================================================================
# 3. PLOT FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3A. Horizontal Z-Score Heatmap  (original layout — unchanged)
#     rows = genes  |  columns = conditions  |  facets = cell types
# ------------------------------------------------------------------------------
generate_zscore_heatmap_horizontal <- function(hm_data, genes, cell_order) {
  
  avg_std <- add_zscore(hm_data) %>%
    group_by(ident, condition, Gene) %>%
    summarise(Mean_Z = mean(Scaled_Expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Gene      = factor(Gene,      levels = rev(genes)),
      ident     = factor(ident,     levels = cell_order),
      condition = factor(condition, levels = c("healed", "not_healed"))
    ) %>%
    filter(!is.na(ident))
  
  lim <- max(abs(avg_std$Mean_Z), na.rm = TRUE)
  
  ggplot(avg_std, aes(x = condition, y = Gene, fill = Mean_Z)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(~ ident, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = CONFIG$zscore_colors["low"], mid = CONFIG$zscore_colors["mid"],
      high = CONFIG$zscore_colors["high"], midpoint = 0,
      limits = c(-lim, lim), name = "Z-Score"
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title              = element_blank(),
      axis.ticks.x.top        = element_blank(),
      axis.ticks.length.x.top = unit(0, "pt"),
      axis.text.x.top         = element_text(angle = 45, hjust = 0, vjust = 0,
                                             size = 9, margin = margin(b = 2)),
      axis.text.y             = element_text(face = "bold.italic", size = 12),
      panel.grid              = element_blank(),
      panel.spacing           = unit(0.1, "lines"),
      strip.placement         = "outside",
      strip.text.x.bottom     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                             face = "bold", size = 9),
      strip.background        = element_blank(),
      legend.position         = "right",
      legend.title            = element_text(size = 11),
      legend.text             = element_text(size = 10)
    )
}


# ------------------------------------------------------------------------------
# 3B. Vertical Z-Score Heatmap  (one narrow panel per gene, side-by-side)
#     rows = cell types  |  columns = healed / not_healed
#
#     Significance asterisks (Wilcoxon rank-sum, healed vs not_healed) are
#     overlaid on the tile for the higher-expressing condition, centred.
#     Empty string = not significant (p >= 0.05).
# ------------------------------------------------------------------------------
generate_zscore_heatmap_vertical <- function(hm_data, genes, cell_order) {
  
  # Compute significance for every (cell type x gene) pair
  sig_df <- compute_wilcoxon(hm_data) %>%
    mutate(ident = factor(ident, levels = rev(cell_order)),
           Gene  = factor(Gene,  levels = genes)) %>%
    filter(!is.na(ident))
  
  # Mean intra-cluster Z-score per (cell type x condition x gene)
  avg_std <- add_zscore(hm_data) %>%
    group_by(ident, condition, Gene) %>%
    summarise(Mean_Z = mean(Scaled_Expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ident     = factor(ident,     levels = rev(cell_order)),
      condition = factor(condition, levels = c("healed", "not_healed")),
      Gene      = factor(Gene,      levels = genes)
    ) %>%
    filter(!is.na(ident))
  
  # Attach significance stars to one representative condition per tile row
  # (stars appear on the "not_healed" column so they are always visible)
  avg_std <- avg_std %>%
    left_join(sig_df %>% select(ident, Gene, stars),
              by = c("ident", "Gene")) %>%
    mutate(
      # Show stars only in the not_healed column to avoid duplication
      label = ifelse(condition == "not_healed", stars, "")
    )
  
  # Shared symmetric colour limit across both genes
  lim <- max(abs(avg_std$Mean_Z), na.rm = TRUE)
  
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  # Build one narrow panel per gene
  panels <- lapply(seq_along(genes), function(i) {
    
    gene_name   <- genes[i]
    df_gene     <- avg_std %>% filter(Gene == gene_name)
    show_legend <- (i == length(genes))
    
    p <- ggplot(df_gene, aes(x = condition, y = ident, fill = Mean_Z)) +
      geom_tile(color = "white", linewidth = 0.6) +
      # Significance asterisks centred in each tile
      geom_text(aes(label = label), size = 4.5, color = "black",
                fontface = "bold", vjust = 0.5) +
      scale_fill_gradient2(
        low = CONFIG$zscore_colors["low"], mid = CONFIG$zscore_colors["mid"],
        high = CONFIG$zscore_colors["high"], midpoint = 0,
        limits = c(-lim, lim), name = "Z-Score"
      ) +
      # Narrow tiles: fix tile width by limiting x expansion
      scale_x_discrete(expand = c(0, 0), labels = cond_labels) +
      scale_y_discrete(expand = c(0, 0)) +
      # Force narrow aspect: each condition column ~0.5 inches wide
      coord_fixed(ratio = 0.5) +
      labs(title = gene_name, x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title      = element_text(face = "bold.italic", size = 14,
                                       hjust = 0.5, margin = margin(b = 6)),
        axis.text.x     = element_text(face = "bold", size = 11,
                                       margin = margin(t = 4)),
        axis.text.y     = element_text(face = "bold", size = 10),
        panel.grid      = element_blank(),
        legend.position = if (show_legend) "right" else "none",
        legend.title    = element_text(size = 11),
        legend.text     = element_text(size = 10)
      )
    
    # Suppress repeated y-axis labels from second panel onward
    if (i > 1) p <- p + theme(axis.text.y = element_blank())
    
    p
  })
  
  # Assemble with shared figure title and subtitle
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = "Intra-cluster Z-Score  \u00b7  healed vs not_healed  \u00b7  * p<0.05  ** p<0.01  *** p<0.001 (Wilcoxon)",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                     margin = margin(b = 8))
      )
    )
}


# ------------------------------------------------------------------------------
# 3C. Custom DotPlot  (one panel per gene, assembled side-by-side)
#
#   Layout per panel:
#     rows  = cell types (ordered top-to-bottom per cell_lineage_order)
#     x     = condition  (healed | not_healed) — plain black labels
#
#   Dot SIZE   = % cells expressing the gene, scaled per gene to [1, 8].
#                Raw % printed as text inside each dot.
#   Dot COLOUR = mean normalised expression (expressors only).
# ------------------------------------------------------------------------------
generate_dotplot_custom <- function(hm_data, genes, cell_order) {
  
  # Aggregate % expressed and mean expression per (cell type x condition x gene)
  dot_data <- hm_data %>%
    group_by(ident, condition, Gene) %>%
    summarise(
      Pct_Expr  = mean(Expression > 0, na.rm = TRUE) * 100,
      Mean_Expr = mean(Expression[Expression > 0], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    replace_na(list(Mean_Expr = 0)) %>%
    mutate(
      ident     = factor(ident,     levels = rev(cell_order)),
      condition = factor(condition, levels = c("healed", "not_healed")),
      Gene      = factor(Gene,      levels = genes)
    ) %>%
    filter(!is.na(ident))
  
  # Per-gene min-max scaling of dot size to [1, 8]
  # Scaling is intentionally per-gene to maximise visual contrast within each
  # gene; raw % is printed inside dots for exact reading.
  dot_data <- dot_data %>%
    group_by(Gene) %>%
    mutate(
      Pct_Scaled = {
        rng <- range(Pct_Expr, na.rm = TRUE)
        if (diff(rng) == 0) rep(4.5, n())
        else 1 + 7 * (Pct_Expr - rng[1]) / diff(rng)
      }
    ) %>%
    ungroup()
  
  # Condition labels — plain text, no colour coding
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  # Build one panel per gene
  panels <- lapply(seq_along(genes), function(i) {
    
    gene_name   <- genes[i]
    df_gene     <- dot_data %>% filter(Gene == gene_name)
    show_legend <- (i == length(genes))
    
    # Build size legend manually using annotated dummy points
    # (avoids the deprecated 'breaks' argument inside guide_legend for
    #  scale_size_identity, which caused warnings in previous version)
    p <- ggplot(df_gene,
                aes(x = condition, y = ident,
                    size  = Pct_Scaled,
                    color = Mean_Expr)) +
      geom_point(alpha = 0.90) +
      scale_color_gradient(
        low  = CONFIG$expr_colors["low"],
        high = CONFIG$expr_colors["high"],
        name = "Mean\nExpression",
        guide = if (show_legend) "colorbar" else "none"
      ) +
      # Use scale_size_identity with a simple guide (no custom breaks)
      scale_size_identity(
        guide = if (show_legend) {
          guide_legend(
            title        = "% Expressed\n(per gene)",
            title.theme  = element_text(size = 10, face = "bold"),
            label.theme  = element_text(size = 9),
            override.aes = list(color = "grey50")
          )
        } else {
          "none"
        }
      ) +
      # Print raw % value inside each dot for exact reading
      geom_text(aes(label = sprintf("%.0f%%", Pct_Expr)),
                size = 2.6, color = "black") +
      # Dashed vertical line separating the two conditions
      geom_vline(xintercept = 1.5, color = "grey75", linewidth = 0.4,
                 linetype = "dashed") +
      scale_x_discrete(labels = cond_labels,
                       expand = expansion(add = 0.65)) +
      scale_y_discrete(expand = expansion(add = 0.55)) +
      labs(title = gene_name, x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title         = element_text(face = "bold.italic", size = 14,
                                          hjust = 0.5, margin = margin(b = 6)),
        # Plain black labels — no colour coding
        axis.text.x        = element_text(face  = "bold", size = 12,
                                          color = "black",
                                          margin = margin(t = 4)),
        axis.text.y        = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = if (show_legend) "right" else "none",
        legend.title       = element_text(size = 10, face = "bold"),
        legend.text        = element_text(size = 9)
      )
    
    # Suppress repeated y-axis labels from second panel onward
    if (i > 1) p <- p + theme(axis.text.y = element_blank())
    
    p
  })
  
  # Assemble panels with shared annotation
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = "Dot size: % expressed (scaled per gene)  \u00b7  Colour: mean normalised expression",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                     margin = margin(b = 8))
      )
    )
}


# ==============================================================================
# 4. PER-PAIR PIPELINE
#    Orchestrates data fetching and all three plot calls for one gene pair.
# ==============================================================================
run_pair_pipeline <- function(seurat_obj, genes, cell_order, output_dir) {
  
  tag <- genes_to_tag(genes)
  message("\n=== Gene pair: ", paste(genes, collapse = " & "), " ===")
  
  # Validate gene names against the Seurat feature space
  valid_genes <- intersect(genes, rownames(seurat_obj))
  missing     <- setdiff(genes, valid_genes)
  if (length(missing) > 0)
    warning("Gene(s) not found in object — skipped: ",
            paste(missing, collapse = ", "))
  if (length(valid_genes) == 0) {
    warning("No valid genes for pair ", tag, ". Aborting.")
    return(invisible(NULL))
  }
  
  # Fetch normalised expression + metadata in one call
  message("  Fetching data...")
  raw_df <- FetchData(seurat_obj,
                      vars  = c(valid_genes, "ident", "condition"),
                      layer = "data")
  
  long_df <- raw_df %>%
    pivot_longer(cols      = all_of(valid_genes),
                 names_to  = "Gene",
                 values_to = "Expression")
  
  # ---- Plot 1: Horizontal Z-Score Heatmap ----------------------------------
  message("  Saving horizontal Z-Score heatmap...")
  p1 <- generate_zscore_heatmap_horizontal(long_df, valid_genes, cell_order)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Horizontal.pdf")),
    plot = p1, width = CONFIG$heatmap_h_width, height = CONFIG$heatmap_h_height,
    dpi = 300
  )
  
  # ---- Plot 2: Vertical Z-Score Heatmap ------------------------------------
  message("  Saving vertical Z-Score heatmap...")
  p2 <- generate_zscore_heatmap_vertical(long_df, valid_genes, cell_order)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Vertical.pdf")),
    plot = p2, width = CONFIG$heatmap_v_width, height = CONFIG$heatmap_v_height,
    dpi = 300
  )
  
  # ---- Plot 3: Custom DotPlot ----------------------------------------------
  message("  Saving custom DotPlot...")
  p3 <- generate_dotplot_custom(long_df, valid_genes, cell_order)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_DotPlot.pdf")),
    plot = p3, width = CONFIG$dotplot_width, height = CONFIG$dotplot_height,
    dpi = 300
  )
  
  message("  Done — 3 PDFs saved  [prefix: ", tag, "]")
  return(invisible(NULL))
}


# ==============================================================================
# 5. MAIN EXECUTION
# ==============================================================================
tryCatch({
  
  message("Setting up output directory...")
  if (!dir.exists(CONFIG$output_dir))
    dir.create(CONFIG$output_dir, recursive = TRUE)
  
  # Load Seurat object — uncomment the line below for actual execution
  # seurat_obj <- readRDS(CONFIG$input_seurat)
  
  # ---- Sanitise metadata ---------------------------------------------------
  message("Sanitising metadata...")
  
  # Normalise condition strings to lowercase
  seurat_obj$condition <- tolower(trimws(seurat_obj$condition))
  
  # Strip numeric cluster prefix from cell-type labels if present
  # (e.g. "3. NK cells" -> "NK cells")
  clean_idents <- str_remove_all(as.character(Idents(seurat_obj)), "^\\d+\\.\\s*")
  Idents(seurat_obj) <- factor(clean_idents, levels = CONFIG$cell_lineage_order)
  
  # Mirror active identity into metadata for FetchData compatibility
  seurat_obj$ident <- Idents(seurat_obj)
  
  # ---- Iterate over gene pairs ---------------------------------------------
  for (pair in CONFIG$gene_pairs) {
    run_pair_pipeline(
      seurat_obj = seurat_obj,
      genes      = pair,
      cell_order = CONFIG$cell_lineage_order,
      output_dir = CONFIG$output_dir
    )
  }
  
  message("\nAll outputs saved to: ", CONFIG$output_dir)
  
}, error = function(e) {
  message("\n[ERROR] Execution halted: ", e$message)
  traceback()
})