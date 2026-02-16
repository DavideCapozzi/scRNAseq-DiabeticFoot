# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)

# 1. Set working directory and file paths
setwd("G:/Drive condivisi/sc-FEDE_DAVIDE/res_validation/violinplots/")

de_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_with_cell_types.xlsx"
seurat_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds"

# ==============================================================================
# PARAMETRI CONFIGURABILI
# ==============================================================================
# Soglia cellule per smoothing aggressivo
min_cells_threshold <- 100 

# Configurazioni specifiche per zoom asse Y (Gene + Cluster)
# Sintassi: "Gene_NomeCluster" = LimiteMassimoY
# Questo serve per zoomare su violini schiacciati senza perdere dati
custom_y_limits <- list(
  "CCL4_CD4 T cells (naive/CM)" = 0.25
)

default_colors <- c("#fdae61", "#7fbc41")
# ==============================================================================

# 2. Load differential expression data
if(!file.exists(de_file)) { 
  stop(paste("ERROR: File not found at:", de_file)) 
}
df <- read_excel(de_file)

# 3. Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(seurat_file)

# 4. Cluster mapping (number -> annotation)
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

# 5. Handle Ident Renaming Robustly
current_idents <- as.character(unique(Idents(seurat_obj)))
if (any(current_idents %in% names(new_names)) || any(current_idents %in% 0:20)) {
  seurat_obj <- RenameIdents(seurat_obj, new_names)
} else {
  current_idents_all <- as.character(seurat_obj@active.ident)
  new_idents <- gsub("^[0-9]+\\. ", "", current_idents_all)
  seurat_obj@active.ident <- factor(new_idents, levels = unique(new_idents))
  names(seurat_obj@active.ident) <- colnames(seurat_obj)
}

# 6. Check and adapt condition column values
cat("\nChecking condition column values...\n")
condition_values <- unique(seurat_obj$condition)
cat("Unique values in 'condition':", paste(condition_values, collapse = ", "), "\n")

if("not_healed" %in% condition_values && "healed" %in% condition_values) {
  my_colors <- c("not_healed" = "#fdae61", "healed" = "#7fbc41")
} else if("Not healed" %in% condition_values && "Healed" %in% condition_values) {
  my_colors <- c("Not healed" = "#fdae61", "Healed" = "#7fbc41")
} else {
  my_colors <- default_colors
  if(length(condition_values) >= 2) names(my_colors) <- condition_values[1:2]
  cat("WARNING: Using default color mapping\n")
}

# 7. Helper function to sanitize filenames
sanitize_filename <- function(name) {
  name <- gsub(" ", "_", name)
  name <- gsub("[\\(\\)/]", "", name)
  name <- gsub("[^A-Za-z0-9_-]", "", name)
  return(name)
}

# 8. Target genes
geni_target <- c("CCL4", "CCR1", "CCR5")

# 9. Counter
total_plots <- 0

# 10. Generate separate violin plot for each gene-cluster combination
for(gene in geni_target) {
  
  cat("\n========================================\n")
  cat("Processing gene:", gene, "\n")
  cat("========================================\n")
  
  if (!gene %in% rownames(seurat_obj)) {
    cat("ERROR: Gene", gene, "not found in Seurat object. Skipping.\n")
    next
  }
  
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
  
  if(length(valid_clusters) == 0) {
    cat("WARNING: No matching clusters found in Seurat object for gene", gene, "\n")
    next
  }
  
  for(cluster in valid_clusters) {
    
    cat("\n  --- Processing cluster:", cluster, "---\n")
    
    # Subset Seurat object
    cells_to_keep <- colnames(seurat_obj)[seurat_obj@active.ident == cluster]
    seurat_subset <- subset(seurat_obj, cells = cells_to_keep)
    
    n_cells <- ncol(seurat_subset)
    cat("  Total cells in cluster:", n_cells, "\n")
    
    # Extract Data directly (ALL DATA, NO FILTERING)
    plot_data <- FetchData(seurat_subset, vars = c(gene, "condition"))
    colnames(plot_data) <- c("expression", "condition")
    
    # Check max expression
    if(max(plot_data$expression) == 0) {
      cat("  WARNING: No expression detected for", gene, "in cluster", cluster, "- Skipping plot.\n")
      next
    }
    
    # Determine smoothing based on cell count
    bandwidth_adjust <- 1.0
    if(n_cells < min_cells_threshold) {
      bandwidth_adjust <- 2.0 
      cat("  Low cell count (<", min_cells_threshold, "). Applying smoothing (adjust=2).\n")
    } else {
      bandwidth_adjust <- 1.2
    }
    
    # CHECK FOR SPECIFIC ZOOM SETTINGS
    # We construct a key "Gene_ClusterName" to check against our custom list
    # If match found, we set a specific Y limit
    current_key <- paste(gene, cluster, sep = "_")
    y_zoom_max <- NULL
    
    if (current_key %in% names(custom_y_limits)) {
      y_zoom_max <- custom_y_limits[[current_key]]
      cat("  APPLYING CUSTOM ZOOM: Y-axis limited to 0 -", y_zoom_max, "for visualization.\n")
    }
    
    safe_cluster_name <- sanitize_filename(cluster)
    pdf_name <- paste0(gene, "_", safe_cluster_name, ".pdf")
    
    # TryCatch handling (Warning handler removed to prevent stopping on minor ggplot warnings)
    tryCatch({
      
      pdf(pdf_name, width = 6, height = 5)
      
      cat("  Generating ggplot violin...\n")
      
      p <- ggplot(plot_data, aes(x = condition, y = expression, fill = condition)) +
        geom_violin(
          scale = "width",      # Forces violin to use max width (Essential for visibility)
          adjust = bandwidth_adjust,
          trim = TRUE,
          linewidth = 0.2       # Correct parameter (replaces 'size')
        ) +
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
          legend.text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank()
        ) +
        ylab("Expression Level")
      
      # APPLY ZOOM IF REQUESTED
      # coord_cartesian zooms the view WITHOUT dropping data points from density calculation
      if (!is.null(y_zoom_max)) {
        p <- p + coord_cartesian(ylim = c(0, y_zoom_max))
      }
      
      print(p)
      dev.off()
      
      cat("  SUCCESS: File generated:", pdf_name, "\n")
      total_plots <- total_plots + 1
      
    }, error = function(e) {
      cat("  ERROR during plotting:", conditionMessage(e), "\n")
      if(dev.cur() != 1) dev.off()
    })
    
  } # End cluster loop
  
} # End gene loop


# ==============================================================================
# NEW SECTION: GLOBAL COMPARISONS & HEATMAP (ROBUST VERSION)
# ==============================================================================

cat("\n========================================\n")
cat("STARTING GLOBAL VISUALIZATION ANALYSIS\n")
cat("========================================\n")

# Create specific output directory
output_dir_global <- "global_comparisons"
if(!dir.exists(output_dir_global)) dir.create(output_dir_global)

# --- 1. PREPARAZIONE DATI COMUNI ---
all_cell_types <- levels(seurat_obj)
n_types <- length(all_cell_types)

# Calcolo larghezza dinamica
# Base 5 pollici + 0.8 pollici per ogni tipo cellulare per spaziatura extra
dyn_width <- 5 + (n_types * 0.8)

# --- 2. GENERAZIONE VIOLIN PLOT ---

for(gene in geni_target) {
  
  cat("Processing global view for:", gene, "\n")
  
  if (!gene %in% rownames(seurat_obj)) {
    cat("  WARNING: Gene", gene, "not found. Skipping.\n")
    next
  }
  
  # Fetch Data
  global_data <- FetchData(seurat_obj, vars = c(gene, "ident", "condition"))
  colnames(global_data) <- c("expression", "cell_type", "condition")
  global_data$cell_type <- factor(global_data$cell_type, levels = levels(seurat_obj))
  
  plot_title <- paste(gene, "\u2212 Expression levels across cell types")
  
  # ----------------------------------------------------------------------------
  # OPTION 1: COMPARABLE SCALES (Shared Y Axis)
  # ----------------------------------------------------------------------------
  pdf_name_comp <- file.path(output_dir_global, paste0(gene, "_AllCells_Comparable_Scale.pdf"))
  
  tryCatch({
    # Height 8 per dare spazio sotto
    pdf(pdf_name_comp, width = dyn_width, height = 8) 
    
    p1 <- ggplot(global_data, aes(x = cell_type, y = expression, fill = condition)) +
      geom_violin(
        scale = "width",       
        position = position_dodge(width = 0.8),
        trim = TRUE,
        linewidth = 0.2
      ) +
      scale_fill_manual(values = my_colors) +
      ggtitle(plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        
        # CRITICAL FIX 1: Increased Left Margin (3cm) and Bottom (3.5cm)
        plot.margin = margin(t = 1, r = 1, b = 3.5, l = 3, "cm")
      ) +
      ylab("Expression Level") +
      coord_cartesian(clip = "off")
    
    print(p1)
    dev.off()
    cat("  -> Generated Comparable Scale plot.\n")
    
  }, error = function(e) {
    cat("  ERROR in Comparable plot:", conditionMessage(e), "\n")
    if(dev.cur() != 1) dev.off()
  })
  
  # ----------------------------------------------------------------------------
  # OPTION 2: MAXIMIZED VISIBILITY (Free Y Scales)
  # ----------------------------------------------------------------------------
  pdf_name_opt <- file.path(output_dir_global, paste0(gene, "_AllCells_Maximized_Shape.pdf"))
  
  tryCatch({
    # Height increased to 9 to allow facets to be taller relative to text
    pdf(pdf_name_opt, width = dyn_width, height = 9)
    
    p2 <- ggplot(global_data, aes(x = condition, y = expression, fill = condition)) +
      geom_violin(
        scale = "width", 
        trim = TRUE, 
        linewidth = 0.2,
        adjust = 1.1
      ) +
      facet_wrap(~cell_type, nrow = 1, scales = "free_y", strip.position = "bottom") + 
      scale_fill_manual(values = my_colors) +
      ggtitle(plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        
        # Strip text (Labels)
        strip.text = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.placement = "outside", # Moves labels outside axis area
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.title.x = element_blank(),
        
        legend.position = "top",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.border = element_blank(),
        axis.line.x = element_line(color="black"),
        
        # CRITICAL FIX 2: Massive Bottom Margin (6cm) to fit long rotated labels
        plot.margin = margin(t = 1, r = 1, b = 6, l = 2, "cm")
      ) +
      coord_cartesian(clip = "off")
    
    print(p2)
    dev.off()
    cat("  -> Generated Maximized Shape plot.\n")
    
  }, error = function(e) {
    cat("  ERROR in Maximized plot:", conditionMessage(e), "\n")
    if(dev.cur() != 1) dev.off()
  })
}

# --- 3. GENERAZIONE HEATMAP TARGET GENES (MANUAL ROBUST METHOD) ---
cat("\nProcessing Target Gene Heatmap...\n")

tryCatch({
  
  # 3a. Manual Aggregation using dplyr
  # This bypasses Seurat's renaming logic entirely.
  
  # Fetch data for genes + metadata
  hm_data <- FetchData(seurat_obj, vars = c(geni_target, "ident", "condition"), slot = "data")
  
  # Rename columns to avoid issues
  # Note: FetchData usually converts dashes to dots in gene names, but here we assume standard names
  # Check if gene names in hm_data match geni_target. If not, map them.
  # Usually Seurat keeps gene names intact in FetchData unless they are non-syntactic.
  
  # Prepare grouping column
  hm_data$Group <- paste(hm_data$ident, hm_data$condition, sep = "_CUTHERE_")
  
  # Helper function for Seurat-style averaging (Log(Mean(Expm1(x)) + 1))
  # Only if working on LogNormalized data (default 'data' slot)
  seurat_mean <- function(x) {
    log1p(mean(expm1(x)))
  }
  
  # Summarize
  library(dplyr)
  avg_df <- hm_data %>%
    group_by(Group) %>%
    summarise(across(all_of(geni_target), seurat_mean))
  
  # 3b. Scale Data (Z-score by Gene)
  # Convert to matrix for scaling
  mat_raw <- as.matrix(avg_df[, -1])
  rownames(mat_raw) <- avg_df$Group
  
  # Scale: We want Z-score PER GENE (Columns in avg_df are genes right now)
  # scale() works on columns, so this is correct for per-gene scaling
  mat_scaled <- scale(mat_raw)
  
  # 3c. Prepare Long Format for ggplot
  # We convert the scaled matrix back to a data frame for plotting
  heatmap_df <- as.data.frame(mat_scaled)
  heatmap_df$Group <- rownames(heatmap_df)
  
  # Pivot Longer manually
  library(tidyr)
  heatmap_long <- pivot_longer(heatmap_df, cols = -Group, names_to = "Gene", values_to = "Expression")
  
  # 3d. Parse Groups back to CellType and Condition
  # We used "_CUTHERE_" as separator to be safe
  split_groups <- strsplit(heatmap_long$Group, "_CUTHERE_")
  heatmap_long$CellType <- sapply(split_groups, `[`, 1)
  heatmap_long$Condition <- sapply(split_groups, `[`, 2)
  
  # 3e. Set Factor Levels for Ordering
  # Gene order
  heatmap_long$Gene <- factor(heatmap_long$Gene, levels = rev(geni_target))
  
  # Cell Type order (Same as Seurat object)
  heatmap_long$CellType <- factor(heatmap_long$CellType, levels = levels(seurat_obj))
  
  # Condition order
  cond_levels <- sort(unique(seurat_obj$condition))
  heatmap_long$Condition <- factor(heatmap_long$Condition, levels = cond_levels)
  
  # 3f. Plotting
  heatmap_filename <- file.path(output_dir_global, "Target_Genes_Heatmap.pdf")
  
  # Calculate width based on number of columns (CellType * Condition)
  n_cols_hm <- length(unique(heatmap_long$Group))
  hm_width <- 5 + (n_cols_hm * 0.4)
  
  pdf(heatmap_filename, width = hm_width, height = 5)
  
  p_hm <- ggplot(heatmap_long, aes(x = interaction(Condition, CellType, sep=" "), y = Gene, fill = Expression)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Z-Score") +
    scale_y_discrete(expand = c(0,0)) +
    # Use nested x-axis logic: Just show CellType on axis, maybe color strip for condition?
    # Keeping it simple: Label is "CellType Condition"
    scale_x_discrete(expand = c(0,0)) +
    ggtitle("Differential Expression of Target Genes") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      # Rotate x labels 90 degrees
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 12, face = "bold.italic"),
      axis.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(t=10, r=10, b=10, l=10)
    )
  
  print(p_hm)
  dev.off()
  
  cat("  -> Generated Target Genes Heatmap (Robust).\n")
  
}, error = function(e) {
  cat("  ERROR in Heatmap generation:", conditionMessage(e), "\n")
  if(dev.cur() != 1) dev.off()
})

cat("\n========================================\n")
cat("Analysis and Visualization Completed.\n")
cat("Outputs saved in:", output_dir_global, "\n")
cat("========================================\n")

cat("\n========================================\n")
cat("Analysis completed successfully\n")
cat("Total plots generated:", total_plots, "\n")
cat("========================================\n")