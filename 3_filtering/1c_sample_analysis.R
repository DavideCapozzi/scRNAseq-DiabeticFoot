
# ==============================================================================
# Script: Conta Cellule per Campione da Oggetto Seurat
# Descrizione: Legge un oggetto Seurat e crea una tabella riassuntiva con il
#              numero di cellule per campione, salvata in formato Excel
# ==============================================================================

# Caricamento librerie necessarie
library(Seurat)
library(writexl)
library(dplyr)

# ==============================================================================
# PARAMETRI DI INPUT/OUTPUT
# ==============================================================================

# Percorso del file RDS contenente l'oggetto Seurat
seurat_obj <- readRDS("G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/2_filtering/merged_seurat_after_QC.rds")

# Directory di output
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/2_filtering/tables"

# Nome del file di output (senza estensione)
output_filename <- "cell_counts_per_sample"

# ==============================================================================
# VERIFICA DIRECTORY DI OUTPUT
# ==============================================================================

# Crea la directory di output se non esiste
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory di output creata:", output_dir, "\n")
}

cat("Oggetto Seurat caricato con successo!\n")
cat("Numero totale di cellule:", ncol(seurat_obj), "\n")

# ==============================================================================
# ESTRAZIONE METADATI E CONTEGGIO CELLULE
# ==============================================================================

# Estrai i metadati
metadata <- seurat_obj@meta.data

# Verifica quale colonna contiene l'informazione sul campione
# Comuni nomi di colonne: "orig.ident", "sample", "Sample", "patient", "donor"
sample_column <- NULL

possible_sample_cols <- c("orig.ident", "sample", "Sample", "patient", 
                          "donor", "sample_id", "SampleID")

for (col in possible_sample_cols) {
  if (col %in% colnames(metadata)) {
    sample_column <- col
    break
  }
}

# Se non trova automaticamente la colonna, usa "orig.ident" come default
if (is.null(sample_column)) {
  sample_column <- "orig.ident"
  cat("Attenzione: usando 'orig.ident' come colonna dei campioni\n")
} else {
  cat("Colonna campioni identificata:", sample_column, "\n")
}

# Conta le cellule per campione
cell_counts <- metadata %>%
  group_by(!!sym(sample_column)) %>%
  summarise(
    N_Cellule = n(),
    .groups = 'drop'
  ) %>%
  rename(Campione = !!sym(sample_column)) %>%
  arrange(Campione)

# Aggiungi riga con il totale
total_row <- data.frame(
  Campione = "TOTALE",
  N_Cellule = sum(cell_counts$N_Cellule)
)

cell_counts_final <- bind_rows(cell_counts, total_row)

# ==============================================================================
# VISUALIZZAZIONE RISULTATI
# ==============================================================================

cat("\n========================================\n")
cat("RIEPILOGO CELLULE PER CAMPIONE\n")
cat("========================================\n\n")
print(cell_counts_final, n = Inf)
cat("\n")

# ==============================================================================
# SALVATAGGIO FILE CSV
# ==============================================================================

csv_path <- paste0(output_dir, "/", output_filename, ".csv")
write.csv(cell_counts_final, 
          file = csv_path, 
          row.names = FALSE,
          quote = FALSE)

cat("File CSV salvato in:", csv_path, "\n")

# ==============================================================================
# CONVERSIONE E SALVATAGGIO IN FORMATO EXCEL
# ==============================================================================

excel_path <- paste0(output_dir, "/", output_filename, ".xlsx")
write_xlsx(cell_counts_final, 
           path = excel_path)

cat("File Excel salvato in:", excel_path, "\n")

# ==============================================================================
# STATISTICHE AGGIUNTIVE (OPZIONALE)
# ==============================================================================

cat("\n========================================\n")
cat("STATISTICHE AGGIUNTIVE\n")
cat("========================================\n")
cat("Numero totale di campioni:", nrow(cell_counts), "\n")
cat("Numero totale di cellule:", sum(cell_counts$N_Cellule), "\n")
cat("Media cellule per campione:", round(mean(cell_counts$N_Cellule), 2), "\n")
cat("Mediana cellule per campione:", median(cell_counts$N_Cellule), "\n")
cat("Min cellule per campione:", min(cell_counts$N_Cellule), "\n")
cat("Max cellule per campione:", max(cell_counts$N_Cellule), "\n")
cat("\n")

cat("Analisi completata con successo!\n")
