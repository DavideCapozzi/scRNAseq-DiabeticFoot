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
#         Tiles annotated with significance asterisks derived from a
#         pre-computed DE results file (FDR = p_val_adj column).
#
#     3.  <TAG>_DotPlot.pdf
#         Custom dot plot — one panel per gene.
#         rows = Seurat clusters (RNA_snn_res.0.7), x = condition.
#         Dot SIZE = % cells expressing the gene (scaled per gene [1,8]).
#         Dot COLOUR = mean normalised expression.
#
#     4.  <TAG>_ViolinPlot_ByCondition.pdf
#         Violin plot — one row of facets per gene (side-by-side), stacked
#         vertically for the two genes.
#         y = normalised expression, x = condition (healed | not_healed).
#         Facets = annotated cell-type labels (same order as heatmaps).
#         Significance asterisks (FDR-based) annotated above each facet.
#
#     5.  <TAG>_ViolinPlot_AllCells.pdf
#         Violin plot — all cells pooled per cell type (healed + not_healed
#         combined), one violin per cell type.
#         y = normalised expression, x = cell type.
#         One panel per gene, stacked vertically.
#
#   Significance thresholds (FDR-based, from DE file — heatmaps & violins):
#     ***  FDR < 0.001
#     **   FDR < 0.01
#     *    FDR < 0.05
#          (no annotation if FDR >= 0.05 or gene not in DE file)
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
  library(readxl)
  library(patchwork)
})


# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------
CONFIG <- list(
  
  # ---- Paths ----------------------------------------------------------------
  input_seurat = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds",
  
  # Pre-computed DE results (healed vs not_healed per cell type).
  # Required columns: gene, p_val_adj, cell_type
  input_de     = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_with_cell_types.xlsx",
  
  output_dir   = "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/heatmaps_dotplots/",
  
  # ---- Gene pairs: one full set of PDFs per pair ---------------------------
  gene_pairs = list(
    c("CCL4", "CCR5"),
    c("CCR1", "CCR8")
  ),
  
  # ---- Colour scales --------------------------------------------------------
  zscore_colors = c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B"),
  expr_colors   = c(low = "#F0F0F0", high = "#B2182B"),
  
  # ---- Biological cell-type order used for heatmaps and violins ------------
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
  heatmap_h_width      = 18,
  heatmap_h_height     = 4,
  heatmap_v_width      = 9,
  heatmap_v_height     = 11,
  dotplot_width        = 14,
  dotplot_height       = 14,
  # Violin by condition: wide to accommodate all cell types side-by-side
  violinplot_cond_width  = 32,
  violinplot_cond_height = 12,
  # Violin all cells: one violin per cell type across x-axis
  violinplot_all_width   = 18,
  violinplot_all_height  = 10
)


# ==============================================================================
# 2. UTILITY FUNCTIONS
# ==============================================================================

# Build a clean file-name prefix from a gene vector (e.g. "CCL4_CCR5")
genes_to_tag <- function(genes) paste(toupper(genes), collapse = "_")

# Add intra-cluster Z-score column ('Scaled_Expr') to a long expression df.
# Scaling is performed independently per (ident, Gene).
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

# Convert a numeric FDR value to a significance asterisk string.
# Returns "***", "**", "*", or "" (not significant / not tested).
fdr_to_stars <- function(fdr) {
  dplyr::case_when(
    is.na(fdr)  ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01  ~ "**",
    fdr < 0.05  ~ "*",
    TRUE        ~ ""
  )
}

# Load and pre-process the DE results file.
# Returns a data frame with columns: gene (uppercase), cell_type, fdr, stars.
# When multiple rows exist for the same (gene, cell_type) pair — e.g. from
# different cluster IDs — the minimum FDR is retained (most conservative
# approach toward false positives while still showing real signals).
load_de_table <- function(path) {
  
  de_raw <- read_excel(path)
  
  required_cols <- c("gene", "p_val_adj", "cell_type")
  missing_cols  <- setdiff(required_cols, colnames(de_raw))
  if (length(missing_cols) > 0)
    stop("DE file is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  
  de_raw %>%
    # Standardise gene names to uppercase for consistent matching
    mutate(
      gene      = toupper(trimws(gene)),
      cell_type = trimws(cell_type),
      p_val_adj = as.numeric(p_val_adj)
    ) %>%
    # Keep the minimum FDR per (gene, cell_type) across all clusters
    group_by(gene, cell_type) %>%
    summarise(fdr = min(p_val_adj, na.rm = TRUE), .groups = "drop") %>%
    mutate(stars = fdr_to_stars(fdr))
}


# ==============================================================================
# 3. PLOT FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3A. Horizontal Z-Score Heatmap  (original layout — unchanged)
#     rows = genes  |  columns = conditions  |  facets = cell types
#     Uses annotated cell-type labels (Idents) and the biological order defined
#     in CONFIG$cell_lineage_order.
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
#     Uses annotated cell-type labels (Idents) and the biological order defined
#     in CONFIG$cell_lineage_order.
#
#     Significance asterisks come from the pre-loaded DE table (FDR-based).
#     Stars are displayed in the centre of the 'not_healed' tile for each
#     (cell type, gene) pair to avoid duplication.
#     Genes absent from the DE file receive no annotation.
# ------------------------------------------------------------------------------
generate_zscore_heatmap_vertical <- function(hm_data, genes, cell_order,
                                             de_table) {
  
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
  
  # Pull significance stars from DE table for the genes in this pair.
  # de_table uses uppercase gene names; match on (gene, cell_type).
  sig_lookup <- de_table %>%
    filter(gene %in% toupper(genes)) %>%
    select(Gene = gene, ident = cell_type, stars)
  
  # Attach stars; default to "" for (gene, cell_type) pairs not in DE file
  avg_std <- avg_std %>%
    mutate(Gene_upper = toupper(as.character(Gene)),
           ident_chr  = as.character(ident)) %>%
    left_join(
      sig_lookup %>% rename(Gene_upper = Gene, ident_chr = ident),
      by = c("Gene_upper", "ident_chr")
    ) %>%
    replace_na(list(stars = "")) %>%
    mutate(
      # Show asterisks only in the not_healed column to avoid duplication
      label = ifelse(condition == "not_healed", stars, "")
    ) %>%
    select(-Gene_upper, -ident_chr)
  
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
      # Significance asterisks centred in the not_healed tile
      geom_text(aes(label = label), size = 4.5, color = "black",
                fontface = "bold", vjust = 0.5) +
      scale_fill_gradient2(
        low = CONFIG$zscore_colors["low"], mid = CONFIG$zscore_colors["mid"],
        high = CONFIG$zscore_colors["high"], midpoint = 0,
        limits = c(-lim, lim), name = "Z-Score"
      ) +
      # coord_fixed keeps tiles narrow (width = 0.5 x height)
      coord_fixed(ratio = 0.5) +
      scale_x_discrete(expand = c(0, 0), labels = cond_labels) +
      scale_y_discrete(expand = c(0, 0)) +
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
  
  # Assemble with shared figure annotation
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = paste0(
        "Intra-cluster Z-Score  \u00b7  healed vs not_healed  ",
        "\u00b7  * FDR<0.05  ** FDR<0.01  *** FDR<0.001"
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                     margin = margin(b = 8))
      )
    )
}


# ------------------------------------------------------------------------------
# 3C. Custom DotPlot — CLUSTER-BASED  (one panel per gene, side-by-side)
#
#   rows  = Seurat clusters (numeric order, top to bottom)
#   x     = condition  (healed | not_healed)
#   Dot SIZE   = % cells expressing the gene, scaled per gene to [1, 8].
#   Dot COLOUR = mean normalised expression (expressors only).
# ------------------------------------------------------------------------------
generate_dotplot_custom <- function(hm_data_clusters, genes, cluster_order) {
  
  # Aggregate % expressed and mean expression per (cluster x condition x gene)
  dot_data <- hm_data_clusters %>%
    group_by(cluster, condition, Gene) %>%
    summarise(
      Pct_Expr  = mean(Expression > 0, na.rm = TRUE) * 100,
      Mean_Expr = mean(Expression[Expression > 0], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    replace_na(list(Mean_Expr = 0)) %>%
    mutate(
      # Use cluster_order for y-axis ordering (reversed so cluster 0 is at top)
      cluster   = factor(cluster,   levels = rev(cluster_order)),
      condition = factor(condition, levels = c("healed", "not_healed")),
      Gene      = factor(Gene,      levels = genes)
    )
  
  # Per-gene min-max scaling of dot size to [1, 8].
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
  
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  # Build one panel per gene
  panels <- lapply(seq_along(genes), function(i) {
    
    gene_name   <- genes[i]
    df_gene     <- dot_data %>% filter(Gene == gene_name)
    show_legend <- (i == length(genes))
    
    p <- ggplot(df_gene,
                aes(x = condition, y = cluster,
                    size  = Pct_Scaled,
                    color = Mean_Expr)) +
      geom_point(alpha = 0.90) +
      scale_color_gradient(
        low   = CONFIG$expr_colors["low"],
        high  = CONFIG$expr_colors["high"],
        name  = "Mean\nExpression",
        guide = if (show_legend) "colorbar" else "none"
      ) +
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
      geom_text(aes(label = sprintf("%.0f%%", Pct_Expr)),
                size = 2.6, color = "black") +
      geom_vline(xintercept = 1.5, color = "grey75", linewidth = 0.4,
                 linetype = "dashed") +
      scale_x_discrete(labels = cond_labels,
                       expand = expansion(add = 0.65)) +
      scale_y_discrete(expand = expansion(add = 0.55)) +
      labs(title = gene_name, x = NULL, y = "Cluster") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title         = element_text(face = "bold.italic", size = 14,
                                          hjust = 0.5, margin = margin(b = 6)),
        axis.text.x        = element_text(face = "bold", size = 12,
                                          color = "black",
                                          margin = margin(t = 4)),
        axis.text.y        = element_text(face = "bold", size = 11),
        axis.title.y       = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = if (show_legend) "right" else "none",
        legend.title       = element_text(size = 10, face = "bold"),
        legend.text        = element_text(size = 9)
      )
    
    if (i > 1) p <- p + theme(axis.text.y = element_blank(),
                              axis.title.y = element_blank())
    p
  })
  
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = "Dot size: % expressed (scaled per gene)  \u00b7  Colour: mean normalised expression  \u00b7  Rows: Seurat clusters (res 0.7)",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                     margin = margin(b = 8))
      )
    )
}


# ------------------------------------------------------------------------------
# 3D. Violin Plot — CELL-TYPE-BASED, BY CONDITION  (publication-ready)
#
#   Publication standards applied (Nature / Cell scRNA-seq style):
#     - Single shared y-axis across all cell-type facets (scales = "fixed").
#       coord_cartesian expands the upper limit to accommodate significance
#       brackets without clipping, while the violin body fills the data range.
#     - All significance brackets drawn at a SINGLE uniform y, so asterisks
#       form a clean horizontal line across the figure.
#     - Fully white background, zero grid lines.
#     - Half-open axis frame (theme_classic): left + bottom lines only.
#     - NO median line, NO boxplot — violin shape alone conveys the distribution.
#     - Strip labels: plain italic text, NO background rectangle.
#     - Condition colours: vivid, colourblind-safe pair.
#     - Legend: horizontal, below the bottom gene panel only.
# ------------------------------------------------------------------------------
generate_violinplot_celltype_condition <- function(hm_data, genes, cell_order,
                                                   de_table) {
  
  # Condition palette — vivid, colourblind-distinguishable
  cond_colors <- c(healed     = "#2196F3",
                   not_healed = "#F44336")
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  # ------------------------------------------------------------------
  # Helper: build significance annotation df for one gene.
  #
  # y_bracket is a SINGLE value computed from the GLOBAL 99th-pctile
  # across ALL cell types × conditions, so every bracket sits at the
  # same height regardless of per-cell-type expression range.
  # y_ceiling is the upper limit passed to coord_cartesian: it adds
  # comfortable room above the bracket text so nothing is clipped.
  # ------------------------------------------------------------------
  build_sig_df <- function(df_gene, gene_name) {
    
    # Hard maximum of the expression data — violin bodies will be capped here.
    # Using max() ensures no violin tip ever crosses the axis line.
    y_max    <- max(df_gene$Expression, na.rm = TRUE)
    y99      <- quantile(df_gene$Expression, 0.99, na.rm = TRUE)
    y_brk    <- y_max * 1.08 + 0.15   # bracket bar: just above the axis cap
    y_txt    <- y_brk + 0.18          # asterisk text height
    y_ceil   <- y_txt + 0.45          # coord_cartesian upper limit (bracket room)
    
    sig_rows <- de_table %>%
      filter(gene == toupper(gene_name), stars != "") %>%
      select(ident = cell_type, stars) %>%
      mutate(ident = factor(ident, levels = cell_order)) %>%
      filter(!is.na(ident))
    
    list(
      sig_df  = if (nrow(sig_rows) == 0) {
        data.frame()
      } else {
        sig_rows %>% mutate(
          y_bracket = y_brk,
          y_text    = y_txt,
          x_lo      = 1, x_hi = 2, x_mid = 1.5
        )
      },
      y_max  = y_max,
      y_ceil = y_ceil
    )
  }
  
  panels <- lapply(seq_along(genes), function(i) {
    
    gene_name   <- genes[i]
    show_legend <- (i == length(genes))
    
    df_gene <- hm_data %>%
      filter(Gene == gene_name) %>%
      mutate(
        ident     = factor(ident,     levels = cell_order),
        condition = factor(condition, levels = c("healed", "not_healed"))
      ) %>%
      filter(!is.na(ident))
    
    sig_out <- build_sig_df(df_gene, gene_name)
    sig_df  <- sig_out$sig_df
    y_max   <- sig_out$y_max
    y_ceil  <- sig_out$y_ceil
    
    p <- ggplot(df_gene,
                aes(x = condition, y = Expression, fill = condition)) +
      
      # ---- Violin body ----
    geom_violin(
      trim      = TRUE,
      scale     = "width",
      alpha     = 1.00,
      linewidth = 0.40,
      color     = "grey20"
    ) +
      
      # ---- Significance bracket at uniform y across all facets ----
    {
      if (nrow(sig_df) > 0) {
        list(
          # Horizontal bar
          geom_segment(
            data        = sig_df,
            aes(x = x_lo, xend = x_hi,
                y = y_bracket, yend = y_bracket),
            inherit.aes = FALSE,
            linewidth   = 0.50,
            color       = "grey10"
          ),
          # Left drop tick
          geom_segment(
            data        = sig_df,
            aes(x = x_lo, xend = x_lo,
                y = y_bracket - 0.06, yend = y_bracket),
            inherit.aes = FALSE,
            linewidth   = 0.50,
            color       = "grey10"
          ),
          # Right drop tick
          geom_segment(
            data        = sig_df,
            aes(x = x_hi, xend = x_hi,
                y = y_bracket - 0.06, yend = y_bracket),
            inherit.aes = FALSE,
            linewidth   = 0.50,
            color       = "grey10"
          ),
          # Asterisk text
          geom_text(
            data        = sig_df,
            aes(x = x_mid, y = y_text, label = stars),
            inherit.aes = FALSE,
            size        = 5,
            fontface    = "bold",
            color       = "grey5",
            vjust       = 0
          )
        )
      }
    } +
      
      # ---- Shared y-axis: hard cap at data max, expand only for brackets ----
    # scale_y_continuous caps violin bodies at the true data maximum so no
    # tip ever crosses the axis line.  coord_cartesian then silently extends
    # the visible area upward to accommodate the significance brackets and
    # asterisk text without clipping (clip = "off").
    facet_wrap(~ ident, nrow = 1, scales = "fixed") +
      scale_y_continuous(limits = c(0, y_max),
                         expand = expansion(mult = c(0, 0))) +
      coord_cartesian(ylim = c(0, y_ceil), clip = "off") +
      
      scale_fill_manual(values = cond_colors, labels = cond_labels,
                        name   = "Condition",
                        guide  = if (show_legend) "legend" else "none") +
      scale_x_discrete(labels = cond_labels) +
      
      labs(title = gene_name, x = NULL, y = "Normalised Expression") +
      
      # ---- Publication theme ----
    theme_classic(base_size = 12) +
      theme(
        plot.title        = element_text(face = "bold.italic", size = 14,
                                         hjust = 0.5, margin = margin(b = 6)),
        axis.text.x       = element_text(angle = 38, hjust = 1, vjust = 1,
                                         size = 8, color = "black"),
        axis.text.y       = element_text(size = 8.5, color = "black"),
        axis.title.y      = element_text(size = 10.5, face = "bold",
                                         margin = margin(r = 5)),
        axis.line         = element_line(color = "black", linewidth = 0.45),
        axis.ticks        = element_line(color = "black", linewidth = 0.38),
        axis.ticks.length = unit(3, "pt"),
        strip.text        = element_text(face = "italic", size = 7.5,
                                         color = "black",
                                         margin = margin(b = 3, t = 2)),
        strip.background  = element_blank(),
        strip.clip        = "off",
        panel.background  = element_rect(fill = "white", color = NA),
        panel.border      = element_blank(),
        panel.grid        = element_blank(),
        panel.spacing     = unit(0.4, "lines"),
        legend.position   = if (show_legend) "bottom" else "none",
        legend.direction  = "horizontal",
        legend.key.size   = unit(9, "pt"),
        legend.title      = element_text(size = 9, face = "bold"),
        legend.text       = element_text(size = 9),
        legend.background = element_blank(),
        legend.margin     = margin(t = 4),
        plot.background   = element_rect(fill = "white", color = NA),
        plot.margin       = margin(4, 6, 4, 6)
      )
    
    p
  })
  
  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = paste0(
        "Normalised expression per annotated cell type  \u00b7  healed vs not_healed  ",
        "\u00b7  * FDR<0.05  ** FDR<0.01  *** FDR<0.001"
      ),
      theme = theme(
        plot.title      = element_text(face = "bold", size = 15, hjust = 0.5,
                                       margin = margin(b = 3)),
        plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey45",
                                       margin = margin(b = 8)),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}


# ------------------------------------------------------------------------------
# 3E. Violin Plot — ALL CELLS POOLED PER CELL TYPE  (publication-ready)
#
#   Layout:
#     - One panel per gene, stacked vertically.
#     - x     = annotated cell type (ordered as CONFIG$cell_lineage_order).
#     - y     = normalised expression (all conditions combined).
#     - Colour = saturated rainbow palette, one hue per cell type.
#               Generated via hcl.colors(n, "Spectral") for perceptual
#               uniformity — vivid, distinct, publication-safe.
#     - NO boxplot overlay; NO grid lines; NO median line or text label.
#     - Half-open axis frame (theme_classic).
#     - Cell-type names on x-axis, rotated 40°, face "bold.italic".
# ------------------------------------------------------------------------------
generate_violinplot_allcells <- function(hm_data, genes, cell_order) {
  
  # Saturated rainbow palette — perceptually uniform via HCL Spectral
  n_ct       <- length(cell_order)
  ct_palette <- setNames(
    hcl.colors(n_ct, palette = "Spectral", rev = FALSE),
    cell_order
  )
  
  panels <- lapply(seq_along(genes), function(i) {
    
    gene_name <- genes[i]
    
    df_gene <- hm_data %>%
      filter(Gene == gene_name) %>%
      mutate(ident = factor(ident, levels = cell_order)) %>%
      filter(!is.na(ident))
    
    ggplot(df_gene, aes(x = ident, y = Expression, fill = ident)) +
      
      # ---- Violin body ----
    geom_violin(
      trim      = TRUE,
      scale     = "width",
      alpha     = 1.00,
      linewidth = 0.40,
      color     = "grey20"
    ) +
      
      scale_fill_manual(values  = ct_palette, guide = "none") +
      scale_color_manual(values = ct_palette, guide = "none") +
      scale_x_discrete(expand = expansion(add = 0.60)) +
      
      labs(title = gene_name, x = NULL, y = "Normalised Expression") +
      
      # ---- Publication theme ----
    theme_classic(base_size = 12) +
      theme(
        plot.title        = element_text(face = "bold.italic", size = 14,
                                         hjust = 0.5, margin = margin(b = 5)),
        axis.text.x       = element_text(angle = 40, hjust = 1, vjust = 1,
                                         size = 9, color = "black",
                                         face = "bold.italic"),
        axis.text.y       = element_text(size = 9, color = "black"),
        axis.title.y      = element_text(size = 10.5, face = "bold",
                                         margin = margin(r = 5)),
        axis.line         = element_line(color = "black", linewidth = 0.45),
        axis.ticks        = element_line(color = "black", linewidth = 0.38),
        axis.ticks.length = unit(3, "pt"),
        panel.background  = element_rect(fill = "white", color = NA),
        panel.border      = element_blank(),
        panel.grid        = element_blank(),
        plot.background   = element_rect(fill = "white", color = NA),
        plot.margin       = margin(4, 6, 4, 6)
      )
  })
  
  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = "Normalised expression per annotated cell type  \u00b7  all cells (healed + not_healed pooled)",
      theme    = theme(
        plot.title      = element_text(face = "bold", size = 15, hjust = 0.5,
                                       margin = margin(b = 3)),
        plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey45",
                                       margin = margin(b = 8)),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}


# ==============================================================================
# 4. PER-PAIR PIPELINE
#    Orchestrates data fetching and all five plot calls for one gene pair.
#
#    Two FetchData calls are performed:
#      - hm_df:      uses annotated cell-type Idents → heatmaps + violins
#      - cluster_df: uses RNA_snn_res.0.7 clusters  → dotplot only
# ==============================================================================
run_pair_pipeline <- function(seurat_obj, genes, cell_order, cluster_order,
                              output_dir, de_table) {
  
  tag <- genes_to_tag(genes)
  message("\n=== Gene pair: ", paste(genes, collapse = " & "), " ===")
  
  # Validate gene names against the Seurat feature space
  valid_genes <- intersect(genes, rownames(seurat_obj))
  missing     <- setdiff(genes, valid_genes)
  if (length(missing) > 0)
    warning("Gene(s) not found in Seurat object — skipped: ",
            paste(missing, collapse = ", "))
  if (length(valid_genes) == 0) {
    warning("No valid genes for pair ", tag, ". Aborting.")
    return(invisible(NULL))
  }
  
  # Warn if some genes are absent from the DE table
  de_genes   <- unique(toupper(de_table$gene))
  missing_de <- setdiff(toupper(valid_genes), de_genes)
  if (length(missing_de) > 0)
    message("  Note: gene(s) not found in DE file (no asterisks): ",
            paste(missing_de, collapse = ", "))
  
  # --------------------------------------------------------------------------
  # Fetch 1: cell-type annotated identity → heatmaps + violin plots
  # --------------------------------------------------------------------------
  message("  Fetching expression data (annotated cell types)...")
  raw_hm <- FetchData(seurat_obj,
                      vars  = c(valid_genes, "ident", "condition"),
                      layer = "data")
  
  hm_df <- raw_hm %>%
    pivot_longer(cols      = all_of(valid_genes),
                 names_to  = "Gene",
                 values_to = "Expression")
  
  # --------------------------------------------------------------------------
  # Fetch 2: cluster identity (RNA_snn_res.0.7) → dotplot
  # --------------------------------------------------------------------------
  message("  Fetching expression data (Seurat clusters — for dotplot)...")
  raw_cl <- FetchData(seurat_obj,
                      vars  = c(valid_genes, "RNA_snn_res.0.7", "condition"),
                      layer = "data")
  
  cluster_df <- raw_cl %>%
    rename(cluster = `RNA_snn_res.0.7`) %>%
    mutate(cluster = as.character(cluster)) %>%
    pivot_longer(cols      = all_of(valid_genes),
                 names_to  = "Gene",
                 values_to = "Expression")
  
  # ---- Plot 1: Horizontal Z-Score Heatmap (cell types) ---------------------
  message("  Saving horizontal Z-Score heatmap...")
  p1 <- generate_zscore_heatmap_horizontal(hm_df, valid_genes, cell_order)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Horizontal.pdf")),
    plot = p1, width = CONFIG$heatmap_h_width, height = CONFIG$heatmap_h_height,
    dpi = 300
  )
  
  # ---- Plot 2: Vertical Z-Score Heatmap (cell types, with FDR asterisks) --
  message("  Saving vertical Z-Score heatmap...")
  p2 <- generate_zscore_heatmap_vertical(hm_df, valid_genes, cell_order,
                                         de_table)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Vertical.pdf")),
    plot = p2, width = CONFIG$heatmap_v_width, height = CONFIG$heatmap_v_height,
    dpi = 300
  )
  
  # ---- Plot 3: Custom DotPlot (clusters) -----------------------------------
  message("  Saving cluster-based DotPlot...")
  p3 <- generate_dotplot_custom(cluster_df, valid_genes, cluster_order)
  ggsave(
    filename = file.path(output_dir, paste0(tag, "_DotPlot.pdf")),
    plot = p3, width = CONFIG$dotplot_width, height = CONFIG$dotplot_height,
    dpi = 300
  )
  
  # ---- Plot 4: Violin by condition (cell types, FDR asterisks) -------------
  message("  Saving cell-type violin plot (healed vs not_healed)...")
  p4 <- generate_violinplot_celltype_condition(hm_df, valid_genes,
                                               cell_order, de_table)
  ggsave(
    filename = file.path(output_dir,
                         paste0(tag, "_ViolinPlot_ByCondition.pdf")),
    plot   = p4,
    width  = CONFIG$violinplot_cond_width,
    height = CONFIG$violinplot_cond_height,
    dpi    = 300
  )
  
  # ---- Plot 5: Violin all cells pooled per cell type -----------------------
  message("  Saving pooled violin plot (all cells per cell type)...")
  p5 <- generate_violinplot_allcells(hm_df, valid_genes, cell_order)
  ggsave(
    filename = file.path(output_dir,
                         paste0(tag, "_ViolinPlot_AllCells.pdf")),
    plot   = p5,
    width  = CONFIG$violinplot_all_width,
    height = CONFIG$violinplot_all_height,
    dpi    = 300
  )
  
  message("  Done — 5 PDFs saved  [prefix: ", tag, "]")
  return(invisible(NULL))
}


# ==============================================================================
# 5. MAIN EXECUTION
# ==============================================================================
tryCatch({
  
  message("Setting up output directory...")
  if (!dir.exists(CONFIG$output_dir))
    dir.create(CONFIG$output_dir, recursive = TRUE)
  
  # ---- Load and validate the DE significance table -------------------------
  message("Loading DE results table...")
  de_table <- load_de_table(CONFIG$input_de)
  message("  DE table loaded: ", nrow(de_table), " unique (gene, cell_type) entries.")
  
  # ---- Load Seurat object --------------------------------------------------
  # Uncomment the line below to load the Seurat object from disk:
  # seurat_obj <- readRDS(CONFIG$input_seurat)
  
  # ---- Sanitise cell-type metadata (used by heatmaps and violins) ----------
  message("Sanitising Seurat cell-type metadata...")
  
  # Normalise condition strings to lowercase
  seurat_obj$condition <- tolower(trimws(seurat_obj$condition))
  
  # Strip numeric cluster prefix from cell-type labels if present
  # (e.g. "3. NK cells" -> "NK cells")
  clean_idents <- str_remove_all(as.character(Idents(seurat_obj)), "^\\d+\\.\\s*")
  Idents(seurat_obj) <- factor(clean_idents, levels = CONFIG$cell_lineage_order)
  
  # Mirror active identity into metadata for FetchData compatibility
  seurat_obj$ident <- Idents(seurat_obj)
  
  # ---- Build cluster order from RNA_snn_res.0.7 ----------------------------
  cluster_order <- sort(
    unique(as.numeric(as.character(seurat_obj$RNA_snn_res.0.7)))
  )
  cluster_order <- as.character(cluster_order)
  message("  Clusters detected (", length(cluster_order), " total): ",
          paste(cluster_order, collapse = ", "))
  
  # ---- Iterate over gene pairs ---------------------------------------------
  for (pair in CONFIG$gene_pairs) {
    run_pair_pipeline(
      seurat_obj    = seurat_obj,
      genes         = pair,
      cell_order    = CONFIG$cell_lineage_order,
      cluster_order = cluster_order,
      output_dir    = CONFIG$output_dir,
      de_table      = de_table
    )
  }
  
  message("\nAll outputs saved to: ", CONFIG$output_dir)
  
}, error = function(e) {
  message("\n[ERROR] Execution halted: ", e$message)
  traceback()
})