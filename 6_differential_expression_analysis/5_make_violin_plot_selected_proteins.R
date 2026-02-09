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
#seurat_obj <- readRDS(seurat_file)

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

cat("\n========================================\n")
cat("Analysis completed successfully\n")
cat("Total plots generated:", total_plots, "\n")
cat("========================================\n")