library(Seurat)
library(dplyr)
library(writexl)

set.seed(123)

# ========================================
# CONFIGURAZIONE DIRECTORY
# ========================================
input_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/2_pre_filtering"   #PATH x CLUSTER FEDE 
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/3_filtering"

# input_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/1_pre_filtering/"
# output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/2_filtering/"

#input_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/2_pre_filtering/"
#output_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/3_filtering/"

# Crea directory output se non esiste
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ========================================
# CARICAMENTO DATI
# ========================================
cat("\n=== LOADING DATA ===\n")
input_file <- file.path(input_dir, "merged_seurat_withAllStats.rds")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

merged_seurat <- readRDS(input_file)
merged_seurat <- JoinLayers(merged_seurat)

cat("Loaded Seurat object from:", input_file, "\n")
cat("Initial cells:", ncol(merged_seurat), "\n")
cat("Initial features:", nrow(merged_seurat), "\n")

# ========================================
# VERIFICA METADATA ESISTENTE
# ========================================
cat("\n=== CHECKING EXISTING QC METRICS ===\n")
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
available_metrics <- qc_metrics[qc_metrics %in% colnames(merged_seurat@meta.data)]
cat("Available QC metrics:", paste(available_metrics, collapse = ", "), "\n")

# Verifica se le statistiche sono già state calcolate
if (!"percent_mito" %in% colnames(merged_seurat@meta.data)) {
  warning("percent_mito not found. This should have been calculated in the previous step!")
  cat("Calculating mitochondrial percentage...\n")
  merged_seurat[["percent_mito"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
}

if (!"percent_ribo" %in% colnames(merged_seurat@meta.data)) {
  warning("percent_ribo not found. This should have been calculated in the previous step!")
  cat("Calculating ribosomal percentage...\n")
  merged_seurat[["percent_ribo"]] <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]")
}

# ========================================
# STATISTICHE PRE-FILTRAGGIO
# ========================================
cat("\n=== PRE-FILTERING STATISTICS ===\n")
cat("Total cells before QC filtering:", ncol(merged_seurat), "\n")

# Statistiche per ogni metrica
cat("\nnCount_RNA (UMI Counts):\n")
cat("  Min:", min(merged_seurat$nCount_RNA), "\n")
cat("  Max:", max(merged_seurat$nCount_RNA), "\n")
cat("  Median:", median(merged_seurat$nCount_RNA), "\n")

cat("\nnFeature_RNA (Genes Detected):\n")
cat("  Min:", min(merged_seurat$nFeature_RNA), "\n")
cat("  Max:", max(merged_seurat$nFeature_RNA), "\n")
cat("  Median:", median(merged_seurat$nFeature_RNA), "\n")

cat("\npercent_mito (Mitochondrial %):\n")
cat("  Min:", round(min(merged_seurat$percent_mito), 2), "%\n")
cat("  Max:", round(max(merged_seurat$percent_mito), 2), "%\n")
cat("  Median:", round(median(merged_seurat$percent_mito), 2), "%\n")

cat("\npercent_ribo (Ribosomal %):\n")
cat("  Min:", round(min(merged_seurat$percent_ribo), 2), "%\n")
cat("  Max:", round(max(merged_seurat$percent_ribo), 2), "%\n")
cat("  Median:", round(median(merged_seurat$percent_ribo), 2), "%\n")

# Conta cellule che non passano ciascun filtro
cells_low_genes <- sum(merged_seurat$nFeature_RNA < 200)
cells_low_count <- sum(merged_seurat$nCount_RNA < 1000)
cells_high_mt <- sum(merged_seurat$percent_mito > 20)
cells_low_ribo <- sum(merged_seurat$percent_ribo < 5)

cat("\nCells failing individual filters:\n")
cat("  nCount_RNA < 1000:", cells_low_count, "(", 
    round(cells_low_count/ncol(merged_seurat)*100, 2), "%)\n")
cat("  nFeature_RNA < 200:", cells_low_genes, "(", 
    round(cells_low_genes/ncol(merged_seurat)*100, 2), "%)\n")
cat("  percent_mito > 20:", cells_high_mt, "(", 
    round(cells_high_mt/ncol(merged_seurat)*100, 2), "%)\n")
cat("  percent_ribo < 5:", cells_low_ribo, "(", 
    round(cells_low_ribo/ncol(merged_seurat)*100, 2), "%)\n")

# ========================================
# APPLICAZIONE FILTRI QC
# ========================================
cat("\n=== APPLYING QC FILTERS ===\n")
cat("Filter criteria:\n")
cat("  - nCount_RNA >= 1000\n")
cat("  - nFeature_RNA >= 200 (no upper limit, doublets already removed)\n")
cat("  - percent_mito <= 20\n")
cat("  - percent_ribo >= 5\n")

# Filtraggio con espressione logica complessa
merged_seurat_filtered <- subset(
  merged_seurat,
  subset = nCount_RNA >= 1000 &
    nFeature_RNA >= 200 & 
    percent_mito <= 20 & 
    percent_ribo >= 5
)

# ========================================
# STATISTICHE POST-FILTRAGGIO
# ========================================
cat("=== POST-FILTERING STATISTICS ===\n")
cells_before <- ncol(merged_seurat)
cells_after <- ncol(merged_seurat_filtered)
cells_removed <- cells_before - cells_after
removal_rate <- (cells_removed / cells_before) * 100

cat("Cells before filtering:", cells_before, "\n")
cat("Cells after filtering:", cells_after, "\n")
cat("Cells removed:", cells_removed, "(", round(removal_rate, 2), "%)\n")
cat("Cells retained:", round(100 - removal_rate, 2), "%\n")

# Statistiche post-filtraggio
cat("\nPost-filtering ranges:\n")
cat("nCount_RNA: [", min(merged_seurat_filtered$nCount_RNA), 
    ", ", max(merged_seurat_filtered$nCount_RNA), "]\n", sep = "")
cat("nFeature_RNA: [", min(merged_seurat_filtered$nFeature_RNA), 
    ", ", max(merged_seurat_filtered$nFeature_RNA), "]\n", sep = "")
cat("percent_mito: [", round(min(merged_seurat_filtered$percent_mito), 2), 
    "%, ", round(max(merged_seurat_filtered$percent_mito), 2), "%]\n", sep = "")
cat("percent_ribo: [", round(min(merged_seurat_filtered$percent_ribo), 2), 
    "%, ", round(max(merged_seurat_filtered$percent_ribo), 2), "%]\n", sep = "")

# ========================================
# VERIFICA INTEGRITÀ DATI
# ========================================
cat("\n=== DATA INTEGRITY CHECK ===\n")
assay_names <- names(Assays(merged_seurat_filtered))
cat("Assays present:", paste(assay_names, collapse = ", "), "\n")

if ("RNA" %in% assay_names) {
  rna_assay <- merged_seurat_filtered[["RNA"]]
  cat("RNA assay slots:", paste(slotNames(rna_assay), collapse = ", "), "\n")
  
  if ("counts" %in% slotNames(rna_assay)) {
    cat("✓ Raw counts preserved\n")
  } else {
    warning("✗ Counts slot missing!")
  }
}

# Verifica metadata essenziali
essential_metadata <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
available <- essential_metadata %in% colnames(merged_seurat_filtered@meta.data)
cat("\nEssential metadata check:\n")
for (i in seq_along(essential_metadata)) {
  status <- ifelse(available[i], "✓", "✗")
  cat(" ", status, essential_metadata[i], "\n")
}

# ========================================
# CARICAMENTO LISTA GENI Y-LINKED E XIST
# ========================================
cat("\n=== LOADING Y-LINKED GENES LIST ===\n")
y_genes_file <- file.path(input_dir, "genes_Y_only_and_xist.txt")

if (!file.exists(y_genes_file)) {
  warning("Y-linked genes file not found: ", y_genes_file)
  warning("Proceeding without Y-chromosome gene filtering")
  y_genes_list <- character(0)
} else {
  y_genes_list <- readLines(y_genes_file)
  y_genes_list <- y_genes_list[nchar(trimws(y_genes_list)) > 0]  # Rimuovi righe vuote
  cat("Total genes in Y-linked list:", length(y_genes_list), "\n")
  cat("Examples from list:", paste(head(y_genes_list, 10), collapse = ", "), "\n")
}

# ========================================
# RIMOZIONE GENI INDESIDERATI
# ========================================
cat("\n=== REMOVING UNWANTED GENES ===\n")

# Salva conteggio iniziale
genes_before <- nrow(merged_seurat_filtered)
cat("Total genes before filtering:", genes_before, "\n\n")

# Pattern per identificare geni da rimuovere (SENZA ribosomiali)
patterns_to_remove <- c(
  "^MT-",           # Mitochondrial genes
  "^HB[^(P|E|S)]",  # Haemoglobin genes (excluding HBP, HBE, HBS)
  "^MALAT1$"        # MALAT1 gene
)

# Identifica geni che matchano i pattern
genes_to_remove <- c()
for (pattern in patterns_to_remove) {
  matched_genes <- grep(pattern, rownames(merged_seurat_filtered), value = TRUE)
  if (length(matched_genes) > 0) {
    cat("Pattern '", pattern, "': ", length(matched_genes), " genes found\n", sep = "")
    cat("  Examples:", paste(head(matched_genes, 5), collapse = ", "), "\n")
    genes_to_remove <- c(genes_to_remove, matched_genes)
  } else {
    cat("Pattern '", pattern, "': 0 genes found\n", sep = "")
  }
}

# Rimuovi duplicati
genes_to_remove <- unique(genes_to_remove)
cat("\nTotal genes matching unwanted patterns:", length(genes_to_remove), "\n")

# ========================================
# IDENTIFICAZIONE E STAMPA GENI Y-LINKED
# ========================================
if (length(y_genes_list) > 0) {
  y_genes_in_data <- intersect(y_genes_list, rownames(merged_seurat_filtered))
  
  cat("\n" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("📋 Y-LINKED GENES FOUND IN DATASET\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  
  if (length(y_genes_in_data) > 0) {
    cat("\nTotal Y-linked genes found:", length(y_genes_in_data), "out of", 
        length(y_genes_list), "in reference list\n")
    cat("Percentage present:", round(length(y_genes_in_data)/length(y_genes_list)*100, 2), "%\n\n")
    
    cat("Complete list of Y-linked genes to be removed:\n")
    cat(paste(rep("-", 80), collapse = ""), "\n", sep = "")
    
    # Stampa in colonne per migliore leggibilità
    genes_per_row <- 5
    for (i in seq(1, length(y_genes_in_data), by = genes_per_row)) {
      end_idx <- min(i + genes_per_row - 1, length(y_genes_in_data))
      cat(sprintf("%-15s", y_genes_in_data[i:end_idx]), "\n")
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n", sep = "")
    
    # Aggiungi ai geni da rimuovere
    genes_to_remove <- c(genes_to_remove, y_genes_in_data)
    
    # Salva lista geni Y trovati in file
    y_genes_output <- file.path(output_dir, "Y_genes_found_in_dataset.txt")
    writeLines(y_genes_in_data, y_genes_output)
    cat("\n✓ List saved to:", y_genes_output, "\n")
    
  } else {
    cat("\n⚠ No Y-linked genes from the reference list found in dataset\n")
  }
  
  cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")
}

# Rimuovi duplicati dopo aggiunta geni Y
genes_to_remove <- unique(genes_to_remove)

# Identifica geni espressi in meno di 5 cellule
gene_counts <- Matrix::rowSums(GetAssayData(merged_seurat_filtered, assay = "RNA", layer = "counts") > 0)
low_expression_genes <- names(gene_counts[gene_counts < 5])
cat("Genes expressed in < 5 cells:", length(low_expression_genes), "\n")

# Combina tutti i geni da rimuovere
all_genes_to_remove <- unique(c(genes_to_remove, low_expression_genes))
cat("\nTotal genes to remove:", length(all_genes_to_remove), "\n")

# Crea breakdown dettagliato per la tabella finale
mt_genes_removed <- length(grep("^MT-", genes_to_remove, value = TRUE))
hb_genes_removed <- length(grep("^HB[^(P|E|S)]", genes_to_remove, value = TRUE))
malat1_removed <- ifelse("MALAT1" %in% genes_to_remove, 1, 0)
y_genes_removed <- ifelse(length(y_genes_list) > 0, 
                          length(intersect(y_genes_list, genes_to_remove)), 
                          0)

# Identifica geni da mantenere
genes_to_keep <- setdiff(rownames(merged_seurat_filtered), all_genes_to_remove)
cat("Genes to retain:", length(genes_to_keep), "\n")

# Applica il filtraggio
merged_seurat_filtered <- subset(merged_seurat_filtered, features = genes_to_keep)

# Statistiche finali
genes_after <- nrow(merged_seurat_filtered)
genes_removed <- genes_before - genes_after
removal_rate_genes <- (genes_removed / genes_before) * 100

cat("\n--- Gene Filtering Summary ---\n")
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Genes removed:", genes_removed, "(", round(removal_rate_genes, 2), "%)\n")
cat("Genes retained:", round(100 - removal_rate_genes, 2), "%\n")

# ========================================
# SALVATAGGIO
# ========================================
cat("\n=== SAVING FILTERED OBJECT ===\n")
output_file <- file.path(output_dir, "merged_seurat_after_QC.rds")
saveRDS(merged_seurat_filtered, output_file)
cat("Saved filtered object to:", output_file, "\n")

# ========================================
# RIEPILOGO FINALE E EXPORT TABELLA
# ========================================
cat("\n" , paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("🎯 QC FILTERING COMPLETED SUCCESSFULLY!\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n", sep = "")

cat("SUMMARY:\n")
cat("✓ QC filters applied on cells:\n")
cat("    • Minimum UMI counts: 1,000\n")
cat("    • Minimum genes: 200\n")
cat("    • Maximum mitochondrial: 20%\n")
cat("    • Minimum ribosomal: 5%\n")
cat("✓ Unwanted genes removed:\n")
cat("    • Mitochondrial genes (^MT-)\n")
cat("    • Haemoglobin genes (^HB[^(P|E|S)])\n")
cat("    • MALAT1 gene\n")
cat("    • Y-linked genes (from custom list)\n")
cat("    • Genes expressed in < 5 cells\n")
cat("    ⚠ Ribosomal genes (^RP[SL]) were NOT removed\n")
cat("✓ Raw counts preserved\n")

# Crea tabella riassuntiva principale
summary_table <- data.frame(
  Metric = c(
    "Total Cells",
    "Cells Removed",
    "Cells Retained (%)",
    "Total Genes",
    "Genes Removed",
    "Genes Retained (%)",
    "Median UMI/Cell",
    "Median Genes/Cell",
    "Median Mito %",
    "Median Ribo %"
  ),
  Before_Filtering = c(
    cells_before,
    "-",
    "100.00",
    genes_before,
    "-",
    "100.00",
    round(median(merged_seurat$nCount_RNA), 0),
    round(median(merged_seurat$nFeature_RNA), 0),
    round(median(merged_seurat$percent_mito), 2),
    round(median(merged_seurat$percent_ribo), 2)
  ),
  After_Filtering = c(
    cells_after,
    cells_removed,
    round(100 - removal_rate, 2),
    genes_after,
    genes_removed,
    round(100 - removal_rate_genes, 2),
    round(median(merged_seurat_filtered$nCount_RNA), 0),
    round(median(merged_seurat_filtered$nFeature_RNA), 0),
    round(median(merged_seurat_filtered$percent_mito), 2),
    round(median(merged_seurat_filtered$percent_ribo), 2)
  ),
  stringsAsFactors = FALSE
)

# Stampa tabella a console
cat("\n=== FILTERING SUMMARY TABLE ===\n")
print(summary_table, row.names = FALSE)

# Crea tabella dettagliata dei geni rimossi per categoria
gene_removal_details <- data.frame(
  Category = c(
    "Mitochondrial (^MT-)",
    "Haemoglobin (^HB[^(P|E|S)])",
    "MALAT1",
    "Y-linked genes",
    "Low expression (< 5 cells)",
    "TOTAL (unique)"
  ),
  Genes_Removed = c(
    mt_genes_removed,
    hb_genes_removed,
    malat1_removed,
    y_genes_removed,
    length(low_expression_genes),
    length(all_genes_to_remove)
  ),
  Percentage_of_Initial = c(
    round(mt_genes_removed / genes_before * 100, 2),
    round(hb_genes_removed / genes_before * 100, 2),
    round(malat1_removed / genes_before * 100, 2),
    round(y_genes_removed / genes_before * 100, 2),
    round(length(low_expression_genes) / genes_before * 100, 2),
    round(length(all_genes_to_remove) / genes_before * 100, 2)
  ),
  stringsAsFactors = FALSE
)

cat("\n=== GENE REMOVAL DETAILS ===\n")
print(gene_removal_details, row.names = FALSE)

# Note sulla rimozione geni ribosomiali
cat("\nNOTE: Ribosomal genes (^RP[SL]) were intentionally KEPT in the dataset.\n")
cat("They are still used for QC filtering (percent_ribo >= 5%) but not removed from the gene list.\n")

# Salva entrambe le tabelle in Excel con fogli multipli
excel_file <- file.path(output_dir, "QC_filtering_summary.xlsx")
write_xlsx(
  list(
    "Summary" = summary_table,
    "Gene_Removal_Details" = gene_removal_details
  ),
  excel_file
)
cat("\n✓ Summary tables saved to:", excel_file, "\n")

cat("\nNEXT STEPS:\n")
cat("  → Normalization (SCTransform or LogNormalize)\n")
cat("  → Scaling and PCA\n")
cat("  → Clustering and UMAP visualization\n")
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
