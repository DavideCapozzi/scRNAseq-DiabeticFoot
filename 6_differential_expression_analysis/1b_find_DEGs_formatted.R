# ==============================================================================
# Script DE - FINAL VERSION (Dual Output)
# ==============================================================================

# --- Libraries ---
library(Seurat)
library(tidyverse)
library(openxlsx)

# ==============================================================================
# 1. SETUP PATHS & DATA
# ==============================================================================

# Input Seurat object path
input_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"

# Output directory
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/6_differential_expression_analysis/results/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# Load Seurat object
message("Loading Seurat object from: ", input_file)
seurat_obj <- readRDS(input_file)

# Set identities
Idents(seurat_obj) <- seurat_obj$RNA_snn_res.0.7
message("Identities set to RNA_snn_res.0.7")

# ==============================================================================
# 2. CORE DE FUNCTION (RETURNS BOTH FORMATS)
# ==============================================================================

run_DE_per_cluster <- function(seurat_obj, cluster_id) {
  
  message(paste("Processing Cluster:", cluster_id))
  
  # 1. Subset cells
  cells_in_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  cluster_subset <- subset(seurat_obj, cells = cells_in_cluster)
  
  # 2. Check conditions
  groups_present <- unique(cluster_subset$outcome)
  
  if(!all(c("healed", "not_healed") %in% groups_present)){
    message(paste("Cluster", cluster_id, "- SKIP: Missing comparison groups. Found:", paste(groups_present, collapse=", ")))
    return(NULL)
  }
  
  # 3. Find Markers (Wilcoxon Test - Single Cell standard)
  de_markers <- FindMarkers(
    cluster_subset,
    ident.1 = "healed",
    ident.2 = "not_healed",
    group.by = "outcome",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    assay = "RNA"
  )
  
  if(nrow(de_markers) == 0){
    message(paste("Cluster", cluster_id, "- No DEGs found."))
    return(NULL)
  }
  
  # Prepare STANDARD format (keep raw output, add Gene and Cluster)
  de_markers_standard <- de_markers %>%
    rownames_to_column("Gene") %>%
    mutate(Cluster = as.character(cluster_id))
  
  # 4. Calculate Average Expression (For Custom Format)
  de_genes <- rownames(de_markers)
  
  avg_exp_obj <- AverageExpression(
    cluster_subset,
    assays = "RNA",
    features = de_genes,
    group.by = "outcome",
    layer = "data", 
    verbose = FALSE
  )
  
  avg_exp_df <- as.data.frame(avg_exp_obj$RNA)
  
  # --- Column Name Handling ---
  cols_found <- colnames(avg_exp_df)
  col_h <- grep("^healed", cols_found, value = TRUE)
  col_nh <- grep("^not.healed", cols_found, value = TRUE) 
  
  if(length(col_h) == 0 || length(col_nh) == 0) {
    warning(paste("Cluster", cluster_id, ": Could not match columns. Found:", paste(cols_found, collapse=", ")))
    return(NULL)
  }
  
  # 5. Build CUSTOM Output Table
  # We start from the standard one but select/rename specific columns
  output_custom <- de_markers_standard %>%
    mutate(
      `Expression Healing` = avg_exp_df[Gene, col_h],
      `Expression Non_healing` = avg_exp_df[Gene, col_nh],
      `Adjusted p-value` = p_val_adj,
      `Significant group` = ifelse(avg_log2FC > 0, "healed", "not_healed")
    ) %>%
    select(
      Gene, 
      `Expression Healing`, 
      `Expression Non_healing`, 
      `Adjusted p-value`, 
      `Significant group`, 
      Cluster
    )
  
  # Return a list containing BOTH formats
  return(list(standard = de_markers_standard, custom = output_custom))
}

# ==============================================================================
# 3. EXECUTION
# ==============================================================================

all_clusters <- levels(Idents(seurat_obj))

if(all(grepl("^[0-9]+$", all_clusters))) {
  all_clusters <- as.character(sort(as.numeric(all_clusters)))
}

message("Starting analysis on ", length(all_clusters), " clusters...")

# Initialize lists to store results
results_custom_list <- list()
results_standard_list <- list()

for(cl in all_clusters){
  res <- run_DE_per_cluster(seurat_obj, cl)
  
  if(!is.null(res)){
    # Separate the results into two lists
    results_standard_list[[cl]] <- res$standard
    results_custom_list[[cl]] <- res$custom
  }
}

# Combine results
final_df_custom <- bind_rows(results_custom_list)
final_df_standard <- bind_rows(results_standard_list)

# ==============================================================================
# 4. SAVING OUTPUTS
# ==============================================================================

# 1. Save CUSTOM format Excel
outfile_custom <- file.path(output_dir, "DE_healed_vs_not_healed_emelieFormat.xlsx")
write.xlsx(final_df_custom, outfile_custom)
message("Custom format results saved to: ", outfile_custom)

# 2. Save STANDARD format Excel (Original format with pct.1, pct.2 etc)
outfile_standard <- file.path(output_dir, "DE_healed_vs_not_healed_ALL_clusts.xlsx")
write.xlsx(final_df_standard, outfile_standard)
message("Standard format results saved to: ", outfile_standard)

# 3. Save SEURAT OBJECT with STANDARD DEGs
# Storing the standard dataframe (with pct info) in misc slot
seurat_obj@misc$DEGs_healed_vs_not_healed <- final_df_standard

seurat_outfile <- file.path(output_dir, "seurat_with_DEGs.rds")
saveRDS(seurat_obj, seurat_outfile)
message("Seurat object (with standard DEGs in @misc) saved to: ", seurat_outfile)