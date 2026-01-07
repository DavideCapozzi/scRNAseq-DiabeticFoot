# ==============================================================================
# Script DE
# ==============================================================================

# --- Libraries ---
library(Seurat)
library(tidyverse)
library(pheatmap)
library(openxlsx)
library(HGNChelper)

# ==============================================================================
# 1. SETUP PATHS & DATA
# ==============================================================================

# Directory dove si trova l'oggetto Seurat
#input_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"
input_file <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/seurat_res_0.7/seurat_res_0.7.rds"

# Directory output
#output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/6_differential_expression_analysis/results"
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/6_differential_expression_analysis/DEGs_res"

# Crea la directory di output se non esiste
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# Caricamento Seurat
message("Loading Seurat object from: ", input_file)
seurat_obj <- readRDS(input_file)

# Imposta identità
Idents(seurat_obj) <- seurat_obj$RNA_snn_res.0.7
message("Identities set to RNA_snn_res.0.7")

# ==============================================================================
# 2. CORE DE FUNCTION
# ==============================================================================

run_DE_per_cluster <- function(seurat_obj, cluster_id) {
  
  # Selezione cellule
  cells_in_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  cluster_subset <- subset(seurat_obj, cells = cells_in_cluster)
  
  # Check presenza gruppi
  if(!all(c("healed", "not_healed") %in% cluster_subset$outcome)){
    message(paste("Cluster", cluster_id, "- SKIP: Missing comparison groups."))
    return(NULL)
  }
  
  # FindMarkers
  de_markers <- FindMarkers(
    cluster_subset,
    ident.1 = "healed",
    ident.2 = "not_healed",
    group.by = "outcome",
    assay = "RNA",
    slot = "data",
    logfc.threshold = 0.25, 
    min.pct = 0.1          
  )
  
  # Formatting
  de_markers <- de_markers %>%
    rownames_to_column("gene") %>%
    mutate(cluster = as.character(cluster_id))
  
  return(de_markers)
}

# ==============================================================================
# 3. RUN ANALYSIS LOOP
# ==============================================================================

all_clusters <- levels(Idents(seurat_obj))
# Sort numerico
if(all(grepl("^[0-9]+$", all_clusters))) {
  all_clusters <- as.character(sort(as.numeric(all_clusters)))
}

message("Starting Differential Expression analysis loop...")
all_markers_list <- lapply(all_clusters, function(cl){
  message("Processing Cluster ", cl, "...")
  run_DE_per_cluster(seurat_obj, cl)
})

# Unione risultati
all_markers_df <- bind_rows(all_markers_list)

# Salvataggio MAIN FILE
outfile_main <- file.path(output_dir, "DE_healed_vs_not_healed_ALL_clusters.xlsx")
write.xlsx(all_markers_df, outfile_main)
message("Main results saved to: ", outfile_main)

# ==============================================================================
# 4. SAVE SEURAT OBJECT WITH DEGs INSIDE
# ==============================================================================

# Aggiungi i DEGs all'oggetto Seurat (slot misc)
seurat_obj@misc$DEGs_healed_vs_not_healed <- all_markers_df

# Definisci path per salvataggio nuovo Seurat
seurat_outfile <- file.path(output_dir, "seurat_with_DEGs.rds")

# Salva l'oggetto
saveRDS(seurat_obj, seurat_outfile)

message("Updated Seurat object saved to: ", seurat_outfile)

