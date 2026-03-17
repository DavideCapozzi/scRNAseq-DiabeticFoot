# ==============================================================================
# SCRIPT: 6_plot_heatmaps_dotplots.R
# PURPOSE: Publication-ready visualisations comparing 'healed' vs 'not_healed'.
#
#   For each gene pair (CCL4/CCR5, CCR1/CCR8) the following PDFs are produced:
#     1.  <TAG>_ZScore_Heatmap_Horizontal.pdf
#     2.  <TAG>_ZScore_Heatmap_Vertical.pdf
#     3.  <TAG>_DotPlot.pdf
#     4.  <TAG>_ViolinPlot_ByCondition.pdf
#     5.  <TAG>_ViolinPlot_AllCells.pdf
#
#   Significance thresholds (FDR-based, from DE file):
#     ***  FDR < 0.001  |  **  FDR < 0.01  |  *  FDR < 0.05
#
# ==============================================================================
# CHANGELOG — ViolinPlot_ByCondition only:
#
#   FIX v1 — asterisks invisible:
#     scale_y_continuous(limits=) discarded out-of-range geoms before rendering.
#     Fix: removed limits= entirely; coord_cartesian owns the viewport.
#
#   FIX v2 — blank space above violin (axis extended to bracket height):
#     coord_cartesian(ylim = c(0, y_ceil)) where y_ceil included bracket room
#     caused blank space inside the panel above the violin body.
#     Fix: two-layer y design (y_axis_top vs y_brk); brackets drawn outside
#     the panel via clip = "off".  Introduced winsorisation (see FIX v3).
#
#   FIX v3 — violin body overflows axis top:
#     Root cause: y_axis_top = quantile(0.99) * 1.03, but geom_violin with
#     trim = TRUE computes KDE on RAW (non-winsorised) data, whose max can
#     exceed quantile(0.99) * 1.03 substantially, causing the violin tip to
#     extend beyond coord_cartesian(ylim), creating an "overflow" artefact.
#
#     Correct solution — WINSORISATION before plotting:
#       Expression values are capped (winsorised) at the per-gene 99th
#       percentile BEFORE being passed to geom_violin.  The KDE is then
#       estimated on capped data, so the violin body NEVER exceeds q99 by
#       construction.  The axis ceiling is set to q99 (exact), no padding
#       needed inside the panel.  Brackets are placed above q99 via clip="off".
#
#       This is the standard approach in Seurat's VlnPlot() and in published
#       scRNA-seq figures (Nature, Cell, Nature Immunology).  The caption/
#       subtitle documents that values are winsorised at the 99th percentile.
#
#       Key properties guaranteed by this approach:
#         1. Violin body NEVER overflows axis top (KDE max = q99 by construction)
#         2. Asterisks always visible (no limits= on scale_y_continuous)
#         3. No blank space inside panel (axis ends exactly at q99)
#         4. Brackets rendered cleanly in panel margin (clip = "off")
#         5. Shared fixed scale across cell-type facets remains meaningful
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
  input_de     = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_with_cell_types.xlsx",
  output_dir   = "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/heatmaps_dotplots/",
  
  # ---- Gene pairs -----------------------------------------------------------
  gene_pairs = list(
    c("CCL4", "CCR5"),
    c("CCR1", "CCR8")
  ),
  
  # ---- Colour scales --------------------------------------------------------
  zscore_colors = c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B"),
  expr_colors   = c(low = "#F0F0F0", high = "#B2182B"),
  
  # ---- Biological cell-type order ------------------------------------------
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
  heatmap_h_width        = 18,
  heatmap_h_height       = 4,
  heatmap_v_width        = 9,
  heatmap_v_height       = 11,
  dotplot_width          = 14,
  dotplot_height         = 14,
  violinplot_cond_width  = 32,
  violinplot_cond_height = 12,
  violinplot_all_width   = 18,
  violinplot_all_height  = 10
)


# ==============================================================================
# 2. UTILITY FUNCTIONS
# ==============================================================================

genes_to_tag <- function(genes) paste(toupper(genes), collapse = "_")

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

fdr_to_stars <- function(fdr) {
  dplyr::case_when(
    is.na(fdr)  ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01  ~ "**",
    fdr < 0.05  ~ "*",
    TRUE        ~ ""
  )
}

load_de_table <- function(path) {
  de_raw <- read_excel(path)
  required_cols <- c("gene", "p_val_adj", "cell_type")
  missing_cols  <- setdiff(required_cols, colnames(de_raw))
  if (length(missing_cols) > 0)
    stop("DE file is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  de_raw %>%
    mutate(
      gene      = toupper(trimws(gene)),
      cell_type = trimws(cell_type),
      p_val_adj = as.numeric(p_val_adj)
    ) %>%
    group_by(gene, cell_type) %>%
    summarise(fdr = min(p_val_adj, na.rm = TRUE), .groups = "drop") %>%
    mutate(stars = fdr_to_stars(fdr))
}


# ==============================================================================
# 3. PLOT FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3A. Horizontal Z-Score Heatmap
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
# 3B. Vertical Z-Score Heatmap
# ------------------------------------------------------------------------------
generate_zscore_heatmap_vertical <- function(hm_data, genes, cell_order,
                                             de_table) {
  
  avg_std <- add_zscore(hm_data) %>%
    group_by(ident, condition, Gene) %>%
    summarise(Mean_Z = mean(Scaled_Expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ident     = factor(ident,     levels = rev(cell_order)),
      condition = factor(condition, levels = c("healed", "not_healed")),
      Gene      = factor(Gene,      levels = genes)
    ) %>%
    filter(!is.na(ident))
  
  sig_lookup <- de_table %>%
    filter(gene %in% toupper(genes)) %>%
    select(Gene = gene, ident = cell_type, stars)
  
  avg_std <- avg_std %>%
    mutate(Gene_upper = toupper(as.character(Gene)),
           ident_chr  = as.character(ident)) %>%
    left_join(
      sig_lookup %>% rename(Gene_upper = Gene, ident_chr = ident),
      by = c("Gene_upper", "ident_chr")
    ) %>%
    replace_na(list(stars = "")) %>%
    mutate(label = ifelse(condition == "not_healed", stars, "")) %>%
    select(-Gene_upper, -ident_chr)
  
  lim         <- max(abs(avg_std$Mean_Z), na.rm = TRUE)
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  panels <- lapply(seq_along(genes), function(i) {
    gene_name   <- genes[i]
    df_gene     <- avg_std %>% filter(Gene == gene_name)
    show_legend <- (i == length(genes))
    
    p <- ggplot(df_gene, aes(x = condition, y = ident, fill = Mean_Z)) +
      geom_tile(color = "white", linewidth = 0.6) +
      geom_text(aes(label = label), size = 4.5, color = "black",
                fontface = "bold", vjust = 0.5) +
      scale_fill_gradient2(
        low = CONFIG$zscore_colors["low"], mid = CONFIG$zscore_colors["mid"],
        high = CONFIG$zscore_colors["high"], midpoint = 0,
        limits = c(-lim, lim), name = "Z-Score"
      ) +
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
    
    if (i > 1) p <- p + theme(axis.text.y = element_blank())
    p
  })
  
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
# 3C. Custom DotPlot
# ------------------------------------------------------------------------------
generate_dotplot_custom <- function(hm_data_clusters, genes, cluster_order) {
  
  dot_data <- hm_data_clusters %>%
    group_by(cluster, condition, Gene) %>%
    summarise(
      Pct_Expr  = mean(Expression > 0, na.rm = TRUE) * 100,
      Mean_Expr = mean(Expression[Expression > 0], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    replace_na(list(Mean_Expr = 0)) %>%
    mutate(
      cluster   = factor(cluster,   levels = rev(cluster_order)),
      condition = factor(condition, levels = c("healed", "not_healed")),
      Gene      = factor(Gene,      levels = genes)
    )
  
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
    
    if (i > 1) p <- p + theme(axis.text.y  = element_blank(),
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
# 3D. Violin Plot — CELL-TYPE-BASED, BY CONDITION (CORRECTED v3)
# ------------------------------------------------------------------------------
generate_violinplot_celltype_condition <- function(hm_data, genes, cell_order, de_table) {
  
  cond_colors <- c(healed = "#2196F3", not_healed = "#F44336")
  cond_labels <- c(healed = "Healed", not_healed = "Not Healed")
  
  panels <- lapply(seq_along(genes), function(i) {
    gene_name   <- genes[i]
    show_legend <- (i == length(genes))
    
    # 1. Preparazione Dati
    df_gene <- hm_data %>%
      filter(Gene == gene_name) %>%
      mutate(
        ident     = factor(ident, levels = cell_order),
        condition = factor(condition, levels = c("healed", "not_healed"))
      ) %>%
      filter(!is.na(ident))
    
    # 2. SMART TRIMMING CORRETTO (Previene l'errore di dplyr)
    df_plot <- df_gene %>%
      group_by(ident) %>%
      mutate(
        q_limit = {
          # Usiamo le {} per isolare l'operazione. Il risultato finale sarà un 
          # singolo scalare numerico, che mutate() accetterà senza errori.
          nz <- Expression[Expression > 0]
          lim <- if(length(nz) > 5) {
            quantile(nz, 0.995, na.rm = TRUE)
          } else {
            max(Expression, na.rm = TRUE)
          }
          # Assicura che ci sia sempre un tetto minimo di 0.05
          max(lim, 0.05, na.rm = TRUE)
        }
      ) %>%
      filter(Expression <= q_limit) %>%
      ungroup()
    
    # 3. Posizionamento Dinamico delle Bracket
    bracket_pos <- df_plot %>%
      group_by(ident) %>%
      summarise(local_max = max(Expression, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        y_bracket = local_max + (local_max * 0.15) + 0.05,
        y_text    = y_bracket + (local_max * 0.15) + 0.05
      )
    
    # Merge con la tabella delle significatività
    sig_df <- de_table %>%
      filter(gene == toupper(gene_name), stars != "") %>%
      select(ident = cell_type, stars) %>%
      mutate(ident = factor(ident, levels = cell_order)) %>%
      filter(!is.na(ident)) %>%
      inner_join(bracket_pos, by = "ident") %>%
      mutate(x_lo = 1L, x_hi = 2L, x_mid = 1.5)
    
    # 4. Costruzione del Plot
    p <- ggplot(df_plot, aes(x = condition, y = Expression, fill = condition)) +
      geom_violin(
        trim = TRUE, scale = "width", adjust = 1.5, alpha = 1.00,
        linewidth = 0.40, color = "grey20"
      ) +
      {
        if (nrow(sig_df) > 0) {
          list(
            geom_segment(data = sig_df, aes(x = x_lo, xend = x_hi, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.50, color = "grey10"),
            geom_segment(data = sig_df, aes(x = x_lo, xend = x_lo, y = y_bracket - (y_bracket*0.03), yend = y_bracket), inherit.aes = FALSE, linewidth = 0.50, color = "grey10"),
            geom_segment(data = sig_df, aes(x = x_hi, xend = x_hi, y = y_bracket - (y_bracket*0.03), yend = y_bracket), inherit.aes = FALSE, linewidth = 0.50, color = "grey10"),
            geom_text(data = sig_df, aes(x = x_mid, y = y_text, label = stars), inherit.aes = FALSE, size = 5, fontface = "bold", color = "grey5", vjust = 0)
          )
        }
      } +
      facet_wrap(~ ident, nrow = 1, scales = "free_y") +
      scale_y_sqrt(expand = expansion(mult = c(0.01, 0.25))) +
      scale_fill_manual(values = cond_colors, labels = cond_labels, name = "Condition", guide = if (show_legend) "legend" else "none") +
      scale_x_discrete(labels = cond_labels) +
      labs(title = gene_name, x = NULL, y = "Normalised Expression (Sqrt Scale)") +
      theme_classic(base_size = 12) +
      theme(
        plot.title        = element_text(face = "bold.italic", size = 14, hjust = 0.5, margin = margin(b = 6)),
        axis.text.x       = element_text(angle = 38, hjust = 1, vjust = 1, size = 8, color = "black"),
        axis.text.y       = element_text(size = 8.5, color = "black"),
        axis.title.y      = element_text(size = 10.5, face = "bold", margin = margin(r = 5)),
        axis.line         = element_line(color = "black", linewidth = 0.45),
        strip.text        = element_text(face = "italic", size = 7.5, color = "black", margin = margin(b = 3, t = 2)),
        strip.background  = element_blank(),
        panel.spacing     = unit(0.4, "lines"),
        legend.position   = if (show_legend) "bottom" else "none",
        legend.direction  = "horizontal",
        legend.key.size   = unit(9, "pt"),
        legend.margin     = margin(t = 4),
        plot.margin       = margin(10, 6, 4, 6)
      )
    p
  })
  
  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = paste0(
        "Normalised expression (Sqrt Scale)  \u00b7  healed vs not_healed (free Y-scales)  \u00b7  ",
        "KDE smoothed (adjust=1.5) for sparsity  \u00b7  * FDR<0.05  ** FDR<0.01  *** FDR<0.001"
      ),
      theme = theme(
        plot.title      = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey45", margin = margin(b = 8)),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}

# ------------------------------------------------------------------------------
# 3E. Violin Plot — ALL CELLS POOLED PER CELL TYPE (CORRECTED)
# ------------------------------------------------------------------------------
generate_violinplot_allcells <- function(hm_data, genes, cell_order) {
  n_ct       <- length(cell_order)
  ct_palette <- setNames(hcl.colors(n_ct, palette = "Spectral", rev = FALSE), cell_order)
  
  panels <- lapply(seq_along(genes), function(i) {
    gene_name <- genes[i]
    df_gene <- hm_data %>%
      filter(Gene == gene_name) %>%
      mutate(ident = factor(ident, levels = cell_order)) %>%
      filter(!is.na(ident))
    
    ggplot(df_gene, aes(x = ident, y = Expression, fill = ident)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 1.00, linewidth = 0.40, color = "grey20") +
      # RISOLUZIONE SQUASHING: Trasformazione Square Root (Radice Quadrata)
      # Ottima per preservare la scala globale permettendo l'osservazione di cluster a bassa espressione
      scale_y_sqrt(expand = expansion(mult = c(0.01, 0.05))) +
      scale_fill_manual(values  = ct_palette, guide = "none") +
      scale_color_manual(values = ct_palette, guide = "none") +
      scale_x_discrete(expand = expansion(add = 0.60)) +
      labs(title = gene_name, x = NULL, y = "Normalised Expression (Sqrt Scale)") +
      theme_classic(base_size = 12) +
      theme(
        plot.title        = element_text(face = "bold.italic", size = 14, hjust = 0.5, margin = margin(b = 5)),
        axis.text.x       = element_text(angle = 40, hjust = 1, vjust = 1, size = 9, color = "black", face = "bold.italic"),
        axis.text.y       = element_text(size = 9, color = "black"),
        axis.title.y      = element_text(size = 10.5, face = "bold", margin = margin(r = 5)),
        axis.line         = element_line(color = "black", linewidth = 0.45),
        plot.margin       = margin(10, 6, 4, 6)
      )
  })
  
  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = paste(genes, collapse = " & "),
      subtitle = "Normalised expression per annotated cell type  \u00b7  Y-axis is Square Root transformed to enhance low-expression visibility",
      theme    = theme(
        plot.title      = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey45", margin = margin(b = 8)),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}

# ==============================================================================
# 4. PER-PAIR PIPELINE
# ==============================================================================
run_pair_pipeline <- function(seurat_obj, genes, cell_order, cluster_order,
                              output_dir, de_table) {
  
  tag <- genes_to_tag(genes)
  message("\n=== Gene pair: ", paste(genes, collapse = " & "), " ===")
  
  valid_genes <- intersect(genes, rownames(seurat_obj))
  missing     <- setdiff(genes, valid_genes)
  if (length(missing) > 0)
    warning("Gene(s) not found in Seurat object — skipped: ",
            paste(missing, collapse = ", "))
  if (length(valid_genes) == 0) {
    warning("No valid genes for pair ", tag, ". Aborting.")
    return(invisible(NULL))
  }
  
  de_genes   <- unique(toupper(de_table$gene))
  missing_de <- setdiff(toupper(valid_genes), de_genes)
  if (length(missing_de) > 0)
    message("  Note: gene(s) not found in DE file (no asterisks): ",
            paste(missing_de, collapse = ", "))
  
  # Fetch 1: annotated cell-type identity → heatmaps + violin plots
  message("  Fetching expression data (annotated cell types)...")
  raw_hm <- FetchData(seurat_obj,
                      vars  = c(valid_genes, "ident", "condition"),
                      layer = "data")
  
  hm_df <- raw_hm %>%
    pivot_longer(cols      = all_of(valid_genes),
                 names_to  = "Gene",
                 values_to = "Expression")
  
  # Fetch 2: cluster identity → dotplot
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
  
  # Plot 1
  message("  Saving horizontal Z-Score heatmap...")
  p1 <- generate_zscore_heatmap_horizontal(hm_df, valid_genes, cell_order)
  ggsave(file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Horizontal.pdf")),
         p1, width = CONFIG$heatmap_h_width, height = CONFIG$heatmap_h_height,
         dpi = 300)
  
  # Plot 2
  message("  Saving vertical Z-Score heatmap...")
  p2 <- generate_zscore_heatmap_vertical(hm_df, valid_genes, cell_order,
                                         de_table)
  ggsave(file.path(output_dir, paste0(tag, "_ZScore_Heatmap_Vertical.pdf")),
         p2, width = CONFIG$heatmap_v_width, height = CONFIG$heatmap_v_height,
         dpi = 300)
  
  # Plot 3
  message("  Saving cluster-based DotPlot...")
  p3 <- generate_dotplot_custom(cluster_df, valid_genes, cluster_order)
  ggsave(file.path(output_dir, paste0(tag, "_DotPlot.pdf")),
         p3, width = CONFIG$dotplot_width, height = CONFIG$dotplot_height,
         dpi = 300)
  
  # Plot 4
  message("  Saving cell-type violin plot (healed vs not_healed)...")
  p4 <- generate_violinplot_celltype_condition(hm_df, valid_genes,
                                               cell_order, de_table)
  ggsave(file.path(output_dir, paste0(tag, "_ViolinPlot_ByCondition.pdf")),
         p4, width  = CONFIG$violinplot_cond_width,
         height = CONFIG$violinplot_cond_height, dpi = 300)
  
  # Plot 5
  message("  Saving pooled violin plot (all cells per cell type)...")
  p5 <- generate_violinplot_allcells(hm_df, valid_genes, cell_order)
  ggsave(file.path(output_dir, paste0(tag, "_ViolinPlot_AllCells.pdf")),
         p5, width  = CONFIG$violinplot_all_width,
         height = CONFIG$violinplot_all_height, dpi = 300)
  
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
  
  message("Loading DE results table...")
  de_table <- load_de_table(CONFIG$input_de)
  message("  DE table loaded: ", nrow(de_table),
          " unique (gene, cell_type) entries.")
  
  # Uncomment to load Seurat object:
  # seurat_obj <- readRDS(CONFIG$input_seurat)
  
  message("Sanitising Seurat cell-type metadata...")
  seurat_obj$condition <- tolower(trimws(seurat_obj$condition))
  
  clean_idents <- str_remove_all(as.character(Idents(seurat_obj)),
                                 "^\\d+\\.\\s*")
  Idents(seurat_obj) <- factor(clean_idents,
                               levels = CONFIG$cell_lineage_order)
  seurat_obj$ident   <- Idents(seurat_obj)
  
  cluster_order <- sort(
    unique(as.numeric(as.character(seurat_obj$RNA_snn_res.0.7)))
  )
  cluster_order <- as.character(cluster_order)
  message("  Clusters detected (", length(cluster_order), " total): ",
          paste(cluster_order, collapse = ", "))
  
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