# LIBRERIE
# ==========================
library(Seurat)
library(ggplot2)
library(dplyr)
library(scds)
library(SingleCellExperiment)
library(scater)

set.seed(123)

# ==========================
# PERCORSI
# ==========================
input_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/0_soupx_analysis/"
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/1_pre_processing/"

# Input Seurat già mergiato
merged_file <- file.path(input_dir, "merged_seurat_soupX_corrected.rds")

# ==========================
# CARICA OGGETTO SEURAT
# ==========================
seurat_obj <- readRDS(merged_file)
cat("Loaded Seurat object with", ncol(seurat_obj), "cells\n\n")

# ==========================
# FUNZIONE PER PROCESSARE UN SINGOLO CAMPIONE
# ==========================
process_sample <- function(seurat_subset, sample_name) {
  initial_cells <- ncol(seurat_subset)
  cat("Processing sample:", sample_name, "\n")
  cat("Cells before doublet removal:", initial_cells, "\n")
  
  # Conversione in SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_subset)
  sce <- logNormCounts(sce)
  
  # Doublet detection con approccio ibrido
  sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE, verb = FALSE)
  
  # Trasferimento risultati a Seurat
  seurat_subset$cxds_score <- sce$cxds_score
  seurat_subset$bcds_score <- sce$bcds_score
  seurat_subset$hybrid_score <- sce$hybrid_score
  seurat_subset$doublet_class <- ifelse(sce$hybrid_call, "doublet", "singlet")
  
  # Statistiche
  n_doublets <- sum(seurat_subset$doublet_class == "doublet")
  doublet_rate <- n_doublets / initial_cells * 100
  expected_rate <- 0.008 * (initial_cells / 1000) * 100

  cat("Doublets detected:", n_doublets, "(", round(doublet_rate, 2), "%)\n")
  cat("Expected doublet rate:", round(expected_rate, 2), "%\n")
  
  # Visualizzazione QC
  print(VlnPlot(seurat_subset, features = "hybrid_score", 
                group.by = "doublet_class", pt.size = 0.1) +
          ggtitle(paste("scds Hybrid -", sample_name, 
                        "\nDoublets:", n_doublets, "(", round(doublet_rate, 2), "%)")) +
          theme_classic())
  
  # Filtraggio doublet
  seurat_subset <- subset(seurat_subset, doublet_class == "singlet")
  filtered_cells <- ncol(seurat_subset)
  cat("Cells after doublet removal:", filtered_cells, "\n\n")
  
  return(seurat_subset)
}

# ==========================
# PROCESSAMENTO PER OGNI CAMPIONE
# ==========================
sample_names <- unique(seurat_obj$sample)
processed_list <- list()

for (sample_name in sample_names) {
  cat("=== Processing sample:", sample_name, "===\n")
  seurat_subset <- subset(seurat_obj, subset = sample == sample_name)
  
  tryCatch({
    processed_list[[sample_name]] <- process_sample(seurat_subset, sample_name)
    cat("✓ Sample processed successfully\n")
  }, error = function(e) {
    cat("✗ ERROR in sample", sample_name, ":", conditionMessage(e), "\n")
  })
  cat("---\n")
}

# ==========================
# UNIONE DI TUTTI I CAMPIONI PROCESSATI
# ==========================
if (length(processed_list) == 1) {
  final_seurat <- processed_list[[1]]
  cat("Only one sample - no merge needed\n")
} else {
  final_seurat <- merge(
    processed_list[[1]],
    y = processed_list[2:length(processed_list)],
    add.cell.ids = names(processed_list)
  )
  cat("Merged", length(processed_list), "samples\n")
}

# ==========================
# AGGIUNTA METADATA PAZIENTI (NUOVA SEZIONE AGGIUNTA)
# ==========================

# Tabella di mapping dai Fastq_File_ID ai nuovi nomi
mapping_table <- data.frame(
  Fastq_File_ID = c("P33554_1001", "P33554_1002", "P33554_1003", "P33554_1004",
                    "P34304_1001", "P34304_1002", "P34304_1003", "P34304_1004",
                    "P34304_1005", "P34304_1007", "P34304_1008", "P34304_1009"),
  outcome = c("not_healed", "not_healed", "healed", "healed",
              "not_healed", "healed", "healed", "healed",
              "not_healed", "not_healed", "healed", "healed"),
  new_name_assigned = c("Patient1", "Patient2", "Patient3", "Patient4",
                       "Patient5", "Patient6", "Patient7", "Patient8",
                       "Patient9", "Patient10", "Patient11", "Patient12")
)

cat("\n=== AGGIUNTA METADATA PAZIENTI ===\n")

# Verifica che tutti i sample siano presenti nella tabella di mapping
current_samples <- unique(final_seurat$sample)
missing_samples <- setdiff(current_samples, mapping_table$Fastq_File_ID)

if (length(missing_samples) > 0) {
  cat("ATTENZIONE: I seguenti sample non sono presenti nella tabella di mapping:\n")
  print(missing_samples)
} else {
  cat("✓ Tutti i sample sono presenti nella tabella di mapping\n")
}

# Aggiungi i metadata all'oggetto Seurat
final_seurat$patient_name <- mapping_table$new_name_assigned[match(final_seurat$sample, mapping_table$Fastq_File_ID)]
final_seurat$outcome <- mapping_table$outcome[match(final_seurat$sample, mapping_table$Fastq_File_ID)]

# Verifica del mapping
cat("\nVerifica del mapping nomi:\n")
verifica_mapping <- final_seurat@meta.data %>%
  select(sample, patient_name, outcome) %>%
  distinct() %>%
  arrange(patient_name)
print(verifica_mapping)

# ==========================
# STATISTICHE FINALI CON NUOVI NOMI
# ==========================
# ==========================
# STATISTICHE FINALI CON NUOVI NOMI
# ==========================
cat("\n")
cat(rep("=", 50), "\n")
cat("STATISTICHE FINALI\n")
cat(rep("=", 50), "\n")

cat("Cellule totali dopo rimozione doublet:", ncol(final_seurat), "\n\n")

cat("Distribuzione per paziente:\n")
print(table(final_seurat$patient_name))
cat("\n")

cat("Distribuzione per outcome:\n")
print(table(final_seurat$outcome))
cat("\n")

# ==========================
# SALVATAGGIO
# ==========================
output_file <- file.path(output_dir, "merged_seurat_after_doublet_QC.rds")
saveRDS(final_seurat, output_file)
cat("Saved processed Seurat object to:", output_file, "\n")

# ==========================
# VERIFICA FINALE
# ==========================
cat("Total cells in final object:", ncol(final_seurat), "\n")
assay_names <- names(Assays(final_seurat))
cat("Assays present:", paste(assay_names, collapse = ", "), "\n")

# Verifica che i nuovi metadata siano stati aggiunti correttamente
cat("\nNuovi metadata aggiunti:\n")
cat("- patient_name:", ifelse("patient_name" %in% names(final_seurat@meta.data), "PRESENTE", "ASSENTE"), "\n")
cat("- outcome:", ifelse("outcome" %in% names(final_seurat@meta.data), "PRESENTE", "ASSENTE"), "\n")

cat("\nProcessamento completato con successo!\n")
