library(Seurat)

seurat_file <- "I:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds"
seurat_obj <- readRDS(seurat_file)

# 1. Costruisci la tabella: conteggio cellule per cluster annotato
cluster_table <- as.data.frame(table(
  Cluster = as.character(Idents(seurat_obj))
))
colnames(cluster_table) <- c("Cluster", "n_cells")
cluster_table$pct <- round(cluster_table$n_cells / sum(cluster_table$n_cells) * 100, 2)
cluster_table <- cluster_table[order(cluster_table$Cluster), ]

# 2. Salva la tabella nel slot misc dell'oggetto (nessun pacchetto extra)
seurat_obj@misc[["cluster_cell_counts"]] <- cluster_table

# 3. Salva l'oggetto aggiornato su disco
output_dir <- "I:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/"  
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(seurat_obj, file = file.path(output_dir, "seurat_obj_with_counts.rds"))

message("Salvato in: ", file.path(output_dir, "seurat_obj_with_counts.rds"))


###########################################

seurat_obj <- readRDS("I:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_obj_with_counts.rds")

# Accesso diretto
seurat_obj@misc[["cluster_cell_counts"]]

# Oppure stampala ordinata per numerosità
seurat_obj@misc[["cluster_cell_counts"]][
  order(-seurat_obj@misc[["cluster_cell_counts"]]$n_cells), 
]