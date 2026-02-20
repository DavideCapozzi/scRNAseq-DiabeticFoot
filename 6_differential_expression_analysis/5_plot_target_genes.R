# ==============================================================================
# SCRIPT: Single Cell Gene Expression Visualization (Violin & Heatmaps)
# AUTHOR: Gemini (Refactored)
# DATE: 2026-02-16
# DESCRIPTION: Generates publication-ready Violin plots (comparable & maximized)
#              and Heatmaps for target genes, comparing conditions.
# ==============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(stringr)

# ==============================================================================
# 1. CENTRALIZED CONFIGURATION
# ==============================================================================
CONFIG <- list(
  # Input Paths
  working_dir = "G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/violinplots_our_sc_data/",
  de_file     = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_with_cell_types.xlsx",
  seurat_file = "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds",
  
  # Analysis Parameters
  genes_target       = c("CCL4", "CCR5"),
  min_cells_smooth   = 100,
  
  # Violin Plot Settings
  cap_outliers       = TRUE, 
  outlier_quantile   = 0.9, 
  
  # General Dimensions (in inches)
  dims = list(
    base_width_violin  = 5,
    width_multiplier   = 0.8,
    height_violin_comp = 8,
    height_violin_max  = 9,
    
    # HEATMAP DIMENSIONS
    base_width_heatmap = 5,      
    height_heatmap_body = 5      # Height of the cells area only (inches)
  ),
  
  # Fonts
  fonts = list(
    title      = 20,
    axis_x     = 10,
    axis_y     = 12,
    strip_text = 11,
    legend     = 12
  ),
  
  # Violin Colors
  colors = c("#fdae61", "#7fbc41"),
  
  # Heatmap Specific Configuration
  heatmap = list(
    # Colors: Blue (Low) -> White (Mid) -> Yellow (High)
    colors = list(
      low  = "#2166ac", 
      mid  = "#FFFFFF", 
      high = "#b2182b"
    ),

    
    # Margins (in cm) - Increased Bottom to 15cm
    margins = list(
      bottom = 12,  
      top    = 1,
      right  = 2,
      left   = 1
    ),
    
    # Font Sizes
    fonts = list(
      condition_label = 9,
      gene_label      = 12,
      cluster_label   = 10
    )
  )
)

# ==============================================================================
# 2. SETUP & DATA LOADING
# ==============================================================================

setwd(CONFIG$working_dir)

# Load Differential Expression Data
if(!file.exists(CONFIG$de_file)) { 
  stop(paste("ERROR: File not found at:", CONFIG$de_file)) 
}
df <- read_excel(CONFIG$de_file)

# Load Seurat Object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(CONFIG$seurat_file)

# Mapping: Cluster IDs to Names
new_names <- c(
  "0" = "T cells (cytotoxic gamma delta)", 
  "1" = "CD4 T cells (CD4 memory/activated Th2 cells)",
  "2" = "CD4 T cells (naive/CM)", 
  "3" = "NK cells (CD16)", 
  "4" = "CD4 T cells (activated CD4 memory T cells)",
  "5" = "CD8 T cells (EM)", 
  "6" = "CD4 T cells (naive/CM)", 
  "7" = "T cells (cytotoxic gamma delta)",
  "8" = "CD4 T cells (naive/CM)", 
  "9" = "B cells (activated/pre-plasmablasts)", 
  "10" = "Classical monocytes", 
  "11" = "B cells (naive/transitional)", 
  "12" = "T cells (NKT-like gamma delta)",
  "13" = "CD4/CD8 T cells (memory T cells)", 
  "14" = "B cells (naive/transitional)", 
  "15" = "Nonclassical monocytes", 
  "16" = "CD4 T cells (naive/CM)", 
  "17" = "Intermediate monocytes", 
  "18" = "pDCs & cDC2"
)

# Robust Renaming Logic
current_idents <- as.character(unique(Idents(seurat_obj)))
if (any(current_idents %in% names(new_names)) || any(current_idents %in% 0:20)) {
  seurat_obj <- RenameIdents(seurat_obj, new_names)
} else {
  # If names already exist but have numeric prefixes (e.g., "0. T cells"), strip them
  current_idents_all <- as.character(seurat_obj@active.ident)
  new_idents <- gsub("^[0-9]+\\. ", "", current_idents_all)
  seurat_obj@active.ident <- factor(new_idents, levels = unique(new_idents))
  names(seurat_obj@active.ident) <- colnames(seurat_obj)
}

# Define Colors based on Condition
cat("\nChecking condition column values...\n")
condition_values <- unique(seurat_obj$condition)
cat("Unique values in 'condition':", paste(condition_values, collapse = ", "), "\n")

if("not_healed" %in% condition_values && "healed" %in% condition_values) {
  my_colors <- c("not_healed" = CONFIG$colors[1], "healed" = CONFIG$colors[2])
} else if("Not healed" %in% condition_values && "Healed" %in% condition_values) {
  my_colors <- c("Not healed" = CONFIG$colors[1], "Healed" = CONFIG$colors[2])
} else {
  my_colors <- CONFIG$colors
  if(length(condition_values) >= 2) names(my_colors) <- condition_values[1:2]
  cat("WARNING: Using default color mapping\n")
}

# Helper: Sanitize Filenames
sanitize_filename <- function(name) {
  name <- gsub(" ", "_", name)
  name <- gsub("[\\(\\)/]", "", name)
  name <- gsub("[^A-Za-z0-9_-]", "", name)
  return(name)
}

# ==============================================================================
# 3. INDIVIDUAL CLUSTER VIOLIN PLOTS (Legacy Logic)
# ==============================================================================
# This section produces one PDF per significant cluster per gene.

total_plots <- 0

# Custom Y limits for specific gene-cluster combinations (Legacy support)
custom_y_limits <- list(
  "CCL4_CD4 T cells (naive/CM)" = 0.25
)

for(gene in CONFIG$genes_target) {
  
  cat("\n========================================\n")
  cat("Processing gene (Individual Plots):", gene, "\n")
  cat("========================================\n")
  
  if (!gene %in% rownames(seurat_obj)) {
    cat("ERROR: Gene", gene, "not found in Seurat object. Skipping.\n")
    next
  }
  
  # Identify significant clusters for this gene
  clusters_sig_for_gene <- df %>%
    filter(gene == !!gene, p_val_adj < 0.05) %>%
    pull(cell_type) %>%
    unique()
  
  if(length(clusters_sig_for_gene) == 0) {
    cat("No significant clusters for", gene, "- skipping\n")
    next
  }
  
  available_clusters <- levels(seurat_obj@active.ident)
  valid_clusters <- intersect(clusters_sig_for_gene, available_clusters)
  
  for(cluster in valid_clusters) {
    
    cat("  --- Processing cluster:", cluster, "---\n")
    
    # Subset Seurat object
    cells_to_keep <- colnames(seurat_obj)[seurat_obj@active.ident == cluster]
    seurat_subset <- subset(seurat_obj, cells = cells_to_keep)
    n_cells <- ncol(seurat_subset)
    
    # Extract Data
    plot_data <- FetchData(seurat_subset, vars = c(gene, "condition"))
    colnames(plot_data) <- c("expression", "condition")
    
    if(max(plot_data$expression) == 0) {
      cat("  WARNING: No expression detected in cluster. Skipping.\n")
      next
    }
    
    # Determine smoothing bandwidth
    bandwidth_adjust <- ifelse(n_cells < CONFIG$min_cells_smooth, 2.0, 1.2)
    
    # Check for custom Y zoom
    current_key <- paste(gene, cluster, sep = "_")
    y_zoom_max <- if (current_key %in% names(custom_y_limits)) custom_y_limits[[current_key]] else NULL
    
    safe_cluster_name <- sanitize_filename(cluster)
    pdf_name <- paste0(gene, "_", safe_cluster_name, ".pdf")
    
    tryCatch({
      pdf(pdf_name, width = 6, height = 5)
      
      p <- ggplot(plot_data, aes(x = condition, y = expression, fill = condition)) +
        geom_violin(scale = "width", adjust = bandwidth_adjust, trim = TRUE, linewidth = 0.2) +
        scale_fill_manual(values = my_colors) +
        ggtitle(paste(gene, "in", cluster)) +
        theme_minimal() + 
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 11, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "top",
          legend.title = element_text(size = 11, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank()
        ) +
        ylab("Expression Level")
      
      if (!is.null(y_zoom_max)) {
        p <- p + coord_cartesian(ylim = c(0, y_zoom_max))
      }
      
      print(p)
      dev.off()
      total_plots <- total_plots + 1
      
    }, error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      if(dev.cur() != 1) dev.off()
    })
  }
}


# ==============================================================================
# 4. GLOBAL COMPARISONS (Violin & Heatmaps)
# ==============================================================================

cat("\n========================================\n")
cat("STARTING GLOBAL VISUALIZATION ANALYSIS\n")
cat("========================================\n")

output_dir_global <- "global_comparisons"
if(!dir.exists(output_dir_global)) dir.create(output_dir_global)

all_cell_types <- levels(seurat_obj)
n_types <- length(all_cell_types)
dyn_width <- CONFIG$dims$base_width_violin + (n_types * CONFIG$dims$width_multiplier)

for(gene in CONFIG$genes_target) {
  
  cat("Processing global view for:", gene, "\n")
  
  if (!gene %in% rownames(seurat_obj)) {
    next
  }
  
  # Fetch Global Data
  global_data <- FetchData(seurat_obj, vars = c(gene, "ident", "condition"))
  colnames(global_data) <- c("expression", "cell_type", "condition")
  global_data$cell_type <- factor(global_data$cell_type, levels = levels(seurat_obj))
  
  plot_title <- paste(gene, "\u2212 Expression levels across cell types")
  
  # ----------------------------------------------------------------------------
  # A. COMPARABLE SCALES (Shared Y Axis)
  # ----------------------------------------------------------------------------
  pdf_name_comp <- file.path(output_dir_global, paste0(gene, "_AllCells_Comparable_Scale.pdf"))
  
  tryCatch({
    pdf(pdf_name_comp, width = dyn_width, height = CONFIG$dims$height_violin_comp) 
    
    p1 <- ggplot(global_data, aes(x = cell_type, y = expression, fill = condition)) +
      geom_violin(scale = "width", position = position_dodge(width = 0.8), trim = TRUE, linewidth = 0.2) +
      scale_fill_manual(values = my_colors) +
      ggtitle(plot_title) +
      theme_minimal() +
      theme(
        plot.title    = element_text(size = CONFIG$fonts$title, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = CONFIG$fonts$axis_x, face = "bold", color = "black"),
        axis.text.y   = element_text(size = CONFIG$fonts$axis_y),
        axis.title.y  = element_text(size = 14, face = "bold"),
        axis.title.x  = element_blank(),
        legend.position = "top",
        legend.text   = element_text(size = CONFIG$fonts$legend),
        panel.grid    = element_blank(),
        axis.line     = element_line(colour = "black"),
        plot.margin   = margin(t = 1, r = 1, b = 3.5, l = 3, "cm")
      ) +
      ylab("Expression Level") +
      coord_cartesian(clip = "off")
    
    print(p1)
    dev.off()
    
  }, error = function(e) {
    cat("  ERROR in Comparable plot:", conditionMessage(e), "\n")
    if(dev.cur() != 1) dev.off()
  })
  
  # ----------------------------------------------------------------------------
  # B. MAXIMIZED VISIBILITY (Free Y Scales + Smart Capping)
  # ----------------------------------------------------------------------------
  pdf_name_opt <- file.path(output_dir_global, paste0(gene, "_AllCells_Maximized_Shape.pdf"))
  
  tryCatch({
    # Prepare data copy for shape maximization
    data_for_shape <- global_data
    
    # OUTLIER HANDLING: Cap values at percentile per group to recover shape
    if (CONFIG$cap_outliers) {
      data_for_shape <- data_for_shape %>%
        group_by(cell_type) %>%
        mutate(
          threshold = quantile(expression[expression > 0], probs = CONFIG$outlier_quantile, na.rm=TRUE),
          expression = ifelse(expression > threshold, threshold, expression)
        ) %>%
        ungroup()
    }
    
    pdf(pdf_name_opt, width = dyn_width, height = CONFIG$dims$height_violin_max)
    
    p2 <- ggplot(data_for_shape, aes(x = condition, y = expression, fill = condition)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, adjust = 1.1) +
      facet_wrap(~cell_type, nrow = 1, scales = "free_y", strip.position = "bottom") + 
      scale_fill_manual(values = my_colors) +
      ggtitle(plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = CONFIG$fonts$title, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        
        # Strip text (Labels)
        strip.text = element_text(angle = 45, hjust = 1, size = CONFIG$fonts$strip_text, face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        
        # Clean axes for "shape only" focus
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.title.x = element_blank(),
        
        legend.position = "top",
        legend.text = element_text(size = CONFIG$fonts$legend),
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        axis.line.x = element_line(color="black"),
        
        # Margin for rotated labels
        plot.margin = margin(t = 1, r = 1, b = 6, l = 2, "cm")
      ) +
      coord_cartesian(clip = "off")
    
    print(p2)
    dev.off()
    
  }, error = function(e) {
    cat("  ERROR in Maximized plot:", conditionMessage(e), "\n")
    if(dev.cur() != 1) dev.off()
  })
}

# ==============================================================================
# 5. HEATMAP (Final Corrected Version - 90 Degree Rotation & Dynamic Folders)
# ==============================================================================
cat("\nProcessing Target Gene Heatmap...\n")

# Clear devices to prevent conflicts
while(!is.null(dev.list())) dev.off()

tryCatch({
  
  # --- 0. Setup Dynamic Output Directory ---
  # Truncate to first 5 genes for folder name to avoid OS path length limits
  genes_for_folder <- if(length(CONFIG$genes_target) > 5) {
    c(CONFIG$genes_target[1:5], "etc")
  } else {
    CONFIG$genes_target
  }
  
  folder_suffix <- paste(genes_for_folder, collapse = "_")
  heatmap_out_dir <- file.path(output_dir_global, paste0("heatmap_", folder_suffix))
  
  if(!dir.exists(heatmap_out_dir)) {
    dir.create(heatmap_out_dir, recursive = TRUE)
  }
  
  cat("  -> Saving heatmaps to:", heatmap_out_dir, "\n")
  
  # --- 1. Fetch Data ---
  hm_data <- FetchData(seurat_obj, vars = c(CONFIG$genes_target, "ident", "condition"), slot = "data")
  
  # --- 2. Prepare Grouping Key ---
  hm_data$Group <- paste(hm_data$ident, hm_data$condition, sep = "_CUTHERE_")
  
  # Helper for log-space averaging
  seurat_mean <- function(x) {
    log1p(mean(expm1(x)))
  }
  
  # --- 3. Aggregation ---
  avg_df <- hm_data %>%
    group_by(Group) %>%
    summarise(across(all_of(CONFIG$genes_target), seurat_mean))
  
  # --- 4. DATA TRANSFORMATION ---
  avg_long <- pivot_longer(avg_df, cols = -Group, names_to = "Gene", values_to = "Expression")
  
  # Parse Metadata
  split_groups <- strsplit(avg_long$Group, "_CUTHERE_")
  avg_long$CellType <- sapply(split_groups, `[`, 1)
  avg_long$Condition <- sapply(split_groups, `[`, 2)
  
  # Set Factors
  avg_long$Gene <- factor(avg_long$Gene, levels = rev(CONFIG$genes_target))
  avg_long$CellType <- factor(avg_long$CellType, levels = levels(seurat_obj))
  avg_long$Condition <- factor(avg_long$Condition, levels = sort(unique(seurat_obj$condition)))
  
  # STRATEGY 1: Standard Global Z-Score
  df_std <- avg_long %>%
    group_by(Gene) %>%
    mutate(Val = as.vector(scale(Expression))) %>%
    ungroup()
  
  # STRATEGY 2: Maximized Relative (Min-Max 0-1)
  df_minmax <- avg_long %>%
    group_by(Gene) %>%
    mutate(Val = (Expression - min(Expression)) / (max(Expression) - min(Expression))) %>%
    ungroup()
  
  # STRATEGY 3: Non-Linear Power Scale
  df_nonlinear <- df_std %>%
    mutate(Val = sign(Val) * sqrt(abs(Val)))
  
  # --- 5. PLOTTING FUNCTION ---
  
  generate_heatmap <- function(data, filename_suffix, fill_scale, plot_title) {
    
    # 1. Width Calculation
    n_celltypes <- length(unique(data$CellType))
    calc_width <- as.numeric(CONFIG$dims$base_width_heatmap)[1] + (n_celltypes * 0.6)
    
    # 2. Height Calculation (Explicit Scalar Fix)
    m_top_inch <- as.numeric(CONFIG$heatmap$margins$top)[1] / 2.54
    m_btm_inch <- as.numeric(CONFIG$heatmap$margins$bottom)[1] / 2.54
    h_body     <- as.numeric(CONFIG$dims$height_heatmap_body)[1]
    
    # Total PDF Height
    calc_height <- h_body + m_top_inch + m_btm_inch
    
    cat(sprintf("  Debug: Plotting %s [W: %.2f, H: %.2f]\n", filename_suffix, calc_width, calc_height))
    
    # Write to the dynamically generated gene-specific directory
    full_filename <- file.path(heatmap_out_dir, paste0("Target_Genes_Heatmap_", filename_suffix, ".pdf"))
    
    pdf(full_filename, width = calc_width, height = calc_height)
    
    p <- ggplot(data, aes(x = Condition, y = Gene, fill = Val)) +
      geom_tile(color = "white") +
      
      # Facet by CellType, label at bottom
      facet_grid(~CellType, switch = "x") + 
      
      fill_scale +
      
      scale_x_discrete(position = "top", expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      
      coord_cartesian(clip = "off") +
      
      ggtitle(plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = as.numeric(CONFIG$fonts$title)[1], face = "bold", hjust = 0.5, margin = margin(b=20)),
        
        # Condition Labels (Top)
        axis.text.x.top = element_text(
          angle = 45, hjust = 0, vjust = 0, 
          size = as.numeric(CONFIG$heatmap$fonts$condition_label)[1], color = "black"
        ),
        axis.title.x = element_blank(),
        
        # Gene Labels (Y Axis)
        axis.text.y = element_text(
          size = as.numeric(CONFIG$heatmap$fonts$gene_label)[1], 
          face = "bold.italic", color = "black"
        ),
        axis.title.y = element_blank(),
        
        # Cluster Labels (Bottom Strip) - ROTATION FIXED TO 90
        strip.text.x.bottom = element_text(
          angle = 90,        
          hjust = 1,         
          vjust = 0.5,       
          size = as.numeric(CONFIG$heatmap$fonts$cluster_label)[1], 
          face = "bold", color = "black"
        ),
        strip.background = element_blank(),
        strip.placement = "outside",
        
        panel.spacing = unit(0.2, "cm"),
        legend.position = "right",
        panel.grid = element_blank(),
        
        # Margins
        plot.margin = margin(
          t = as.numeric(CONFIG$heatmap$margins$top)[1], 
          r = as.numeric(CONFIG$heatmap$margins$right)[1], 
          b = as.numeric(CONFIG$heatmap$margins$bottom)[1], 
          l = as.numeric(CONFIG$heatmap$margins$left)[1], 
          unit = "cm"
        )
      )
    
    print(p)
    dev.off()
    cat("  -> Generated:", filename_suffix, "\n")
  }
  
  # --- 6. EXECUTE PLOTS ---
  
  col_low  <- CONFIG$heatmap$colors$low
  col_mid  <- CONFIG$heatmap$colors$mid
  col_high <- CONFIG$heatmap$colors$high
  
  # Calculate dynamic symmetric limits for divergent scales to balance the legend
  max_z <- max(abs(df_std$Val), na.rm = TRUE)
  if(max_z == 0 || !is.finite(max_z)) max_z <- 1 # Fallback for edge cases
  
  max_pwr <- max(abs(df_nonlinear$Val), na.rm = TRUE)
  if(max_pwr == 0 || !is.finite(max_pwr)) max_pwr <- 1
  
  # 1. Standard Z-Score
  generate_heatmap(
    df_std, "1_Standard_ZScore",
    scale_fill_gradient2(
      low = col_low, mid = col_mid, high = col_high, midpoint = 0, 
      limits = c(-max_z, max_z), name = "Z-Score"
    ),
    "Differential Expression (Standard Z-Score)"
  )
  
  # 2. Maximized Relative
  generate_heatmap(
    df_minmax, "2_Maximized_Relative",
    scale_fill_gradientn(
      colors = c(col_low, col_mid, col_high), values = c(0, 0.5, 1), 
      limits = c(0, 1), name = "Rel. Expr\n(0-1)"
    ),
    "Differential Expression (Relative Min-Max)"
  )
  
  # 3. Non-Linear Power Scale
  generate_heatmap(
    df_nonlinear, "3_NonLinear_PowerScale",
    scale_fill_gradient2(
      low = col_low, mid = col_mid, high = col_high, midpoint = 0, 
      limits = c(-max_pwr, max_pwr), name = "Power(Z)"
    ),
    "Differential Expression (Gamma Corrected)"
  )
  
}, error = function(e) {
  cat("  ERROR in Heatmap generation:", conditionMessage(e), "\n")
  if(dev.cur() != 1) dev.off()
})

cat("\n========================================\n")
cat("Analysis completed successfully\n")
cat("Outputs saved in:", output_dir_global, "and its subdirectories\n")
cat("Total individual plots generated:", total_plots, "\n")
cat("========================================\n")