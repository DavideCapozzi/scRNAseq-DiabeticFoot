library(Seurat)
library(dplyr)
library(writexl)

set.seed(123)

# ========================================
# CONFIGURAZIONE DIRECTORY
# ========================================
input_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/2_pre_filtering"
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/3_filtering"

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

if (!"percent_mito" %in% colnames(merged_seurat@meta.data)) {
  warning("percent_mito not found!")
  merged_seurat[["percent_mito"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
}

if (!"percent_ribo" %in% colnames(merged_seurat@meta.data)) {
  warning("percent_ribo not found!")
  merged_seurat[["percent_ribo"]] <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]")
}

# ========================================
# STATISTICHE PRE-FILTRAGGIO
# ========================================
cat("\n=== PRE-FILTERING STATISTICS ===\n")
cat("Total cells before QC filtering:", ncol(merged_seurat), "\n")

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
cat("  - nFeature_RNA >= 200\n")
cat("  - percent_mito <= 20\n")
cat("  - percent_ribo >= 5\n")

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
cat("\n=== POST-FILTERING STATISTICS ===\n")
cells_before <- ncol(merged_seurat)
cells_after <- ncol(merged_seurat_filtered)
cells_removed <- cells_before - cells_after
removal_rate <- (cells_removed / cells_before) * 100

cat("Cells before filtering:", cells_before, "\n")
cat("Cells after filtering:", cells_after, "\n")
cat("Cells removed:", cells_removed, "(", round(removal_rate, 2), "%)\n")
cat("Cells retained:", round(100 - removal_rate, 2), "%\n")

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
  cat("Layer names in RNA assay:", paste(names(rna_assay@layers), collapse = ", "), "\n")
}

essential_metadata <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
available <- essential_metadata %in% colnames(merged_seurat_filtered@meta.data)
cat("\nEssential metadata check:\n")
for (i in seq_along(essential_metadata)) {
  status <- ifelse(available[i], "✓", "✗")
  cat(" ", status, essential_metadata[i], "\n")
}

# ========================================
# CARICAMENTO LISTA GENI Y-LINKED
# ========================================
cat("\n=== LOADING Y-LINKED GENES LIST ===\n")
y_genes_file <- file.path(input_dir, "genes_Y_only_and_xist.txt")

if (!file.exists(y_genes_file)) {
  warning("Y-linked genes file not found")
  y_genes_list <- character(0)
} else {
  y_genes_list <- readLines(y_genes_file)
  y_genes_list <- y_genes_list[nchar(trimws(y_genes_list)) > 0]
  cat("Total genes in Y-linked list:", length(y_genes_list), "\n")
  cat("Examples from list:", paste(head(y_genes_list, 10), collapse = ", "), "\n")
}

# ========================================
# RIMOZIONE GENI INDESIDERATI
# ========================================
cat("\n=== REMOVING UNWANTED GENES ===\n")

genes_before <- nrow(merged_seurat_filtered)
cat("Total genes before filtering:", genes_before, "\n\n")

patterns_to_remove <- c(
  "^MT-",
  "^HB[^(P|E|S)]",
  "^MALAT1$"
)

genes_to_remove <- c()
for (pattern in patterns_to_remove) {
  matched_genes <- grep(pattern, rownames(merged_seurat_filtered), value = TRUE)
  if (length(matched_genes) > 0) {
    cat("Pattern '", pattern, "': ", length(matched_genes), " genes found\n", sep = "")
    cat("  Examples:", paste(head(matched_genes, 5), collapse = ", "), "\n")
    genes_to_remove <- c(genes_to_remove, matched_genes)
  }
}

genes_to_remove <- unique(genes_to_remove)
cat("\nTotal genes matching unwanted patterns:", length(genes_to_remove), "\n")

# ========================================
# IDENTIFICAZIONE GENI Y-LINKED
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
    
    genes_per_row <- 5
    for (i in seq(1, length(y_genes_in_data), by = genes_per_row)) {
      end_idx <- min(i + genes_per_row - 1, length(y_genes_in_data))
      cat(sprintf("%-15s", y_genes_in_data[i:end_idx]), "\n")
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n", sep = "")
    genes_to_remove <- c(genes_to_remove, y_genes_in_data)
    
    y_genes_output <- file.path(output_dir, "Y_genes_found_in_dataset.txt")
    writeLines(y_genes_in_data, y_genes_output)
    cat("\n✓ List saved to:", y_genes_output, "\n")
  } else {
    cat("\n⚠ No Y-linked genes found in dataset\n")
  }
  
  cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")
}

genes_to_remove <- unique(genes_to_remove)

# ========================================
# IDENTIFICARE GENE COUNTS
# ========================================
cat("\n=== IDENTIFYING GENE COUNTS (Multiple Layer Seurat v5) ===\n")

rna_assay <- merged_seurat_filtered[["RNA"]]
layer_names <- names(rna_assay@layers)
cat("Number of layers:", length(layer_names), "\n")
cat("Layer names:", paste(layer_names, collapse = ", "), "\n\n")

gene_counts <- Matrix::rowSums(rna_assay@layers[[1]] > 0)

for (i in 2:length(layer_names)) {
  gene_counts <- gene_counts + Matrix::rowSums(rna_assay@layers[[i]] > 0)
}

cat("Gene counts calculated from", length(layer_names), "layers\n")
cat("Total cells with expression data:", sum(gene_counts), "\n")

low_expression_genes <- names(gene_counts[gene_counts < 5])
cat("Genes expressed in < 5 cells:", length(low_expression_genes), "\n")

all_genes_to_remove <- unique(c(genes_to_remove, low_expression_genes))
cat("\nTotal genes to remove:", length(all_genes_to_remove), "\n")

mt_genes_removed <- length(grep("^MT-", genes_to_remove, value = TRUE))
hb_genes_removed <- length(grep("^HB[^(P|E|S)]", genes_to_remove, value = TRUE))
malat1_removed <- ifelse("MALAT1" %in% genes_to_remove, 1, 0)
y_genes_removed <- ifelse(length(y_genes_list) > 0, 
                          length(intersect(y_genes_list, genes_to_remove)), 
                          0)

genes_to_keep <- setdiff(rownames(merged_seurat_filtered), all_genes_to_remove)
cat("Genes to retain:", length(genes_to_keep), "\n")

merged_seurat_filtered <- subset(merged_seurat_filtered, features = genes_to_keep)

genes_after <- nrow(merged_seurat_filtered)
genes_removed <- genes_before - genes_after
removal_rate_genes <- (genes_removed / genes_before) * 100

cat("\n--- Gene Filtering Summary ---\n")
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Genes removed:", genes_removed, "(", round(removal_rate_genes, 2), "%)\n")
cat("Genes retained:", round(100 - removal_rate_genes, 2), "%\n")

# ========================================
# UNIFICAZIONE LAYERS
# ========================================
cat("\n=== UNIFYING LAYERS INTO SINGLE COUNT MATRIX ===\n")
cat("Current number of layers:", length(names(merged_seurat_filtered[["RNA"]]@layers)), "\n")

rna_assay <- merged_seurat_filtered[["RNA"]]
layer_names <- names(rna_assay@layers)

# Concatena orizzontalmente i layer (cbind) perché hanno dimensioni diverse
cat("Layer dimensions:\n")
for (i in 1:length(layer_names)) {
  cat("  ", layer_names[i], ": ", dim(rna_assay@layers[[i]])[1], " x ", 
      dim(rna_assay@layers[[i]])[2], " genes x cells\n", sep = "")
}

# Concatenazione orizzontale (combina le cellule)
combined_counts <- rna_assay@layers[[1]]
for (i in 2:length(layer_names)) {
  combined_counts <- cbind(combined_counts, rna_assay@layers[[i]])
}

cat("\nCombined count matrix dimensions:", dim(combined_counts), "\n")
cat("✓ Successfully combined", length(layer_names), "layers\n")

merged_seurat_filtered[["RNA"]] <- CreateAssay5Object(counts = combined_counts)

cat("✓ Layers unified into single 'counts' layer\n")
cat("New number of layers:", length(names(merged_seurat_filtered[["RNA"]]@layers)), "\n")

# ========================================
# RIMOZIONE METADATA NON UTILI / RIDENOMINAZIONE
# ========================================
cat("\n=== CLEANING UP METADATA ===\n")

# Mantieni solo i metadata importanti per le analisi
metadata_cols_to_keep <- c(
  "orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo",
  "sample", "patient_name", "outcome", "condition"
)

# Verifica quali colonne esistono
metadata_cols_available <- metadata_cols_to_keep[metadata_cols_to_keep %in% colnames(merged_seurat_filtered@meta.data)]

# Mantieni solo quelle disponibili
merged_seurat_filtered@meta.data <- merged_seurat_filtered@meta.data[, metadata_cols_available, drop = FALSE]

cat("Kept metadata columns:", paste(metadata_cols_available, collapse = ", "), "\n")
cat("Total metadata columns:", ncol(merged_seurat_filtered@meta.data), "\n")

# ========================================
# VERIFICA FINALE PRE-SALVATAGGIO
# ========================================
cat("\n=== FINAL VERIFICATION ===\n")
cat("Final object dimensions: ", nrow(merged_seurat_filtered), " genes x ", 
    ncol(merged_seurat_filtered), " cells\n", sep = "")
cat("Memory usage (approximate):", format(object.size(merged_seurat_filtered), units = "MB"), "\n")

# Controlla che non ci siano NA o valori strani
cat("Checking for NA values in counts matrix...\n")
if (any(is.na(combined_counts@x))) {
  warning("NA values found in counts matrix!")
} else {
  cat("✓ No NA values in counts\n")
}

# Verifica che tutti i samples siano ancora presenti
cat("Samples present:", paste(unique(merged_seurat_filtered$orig.ident), collapse = ", "), "\n")

# ========================================
# SALVATAGGIO
# ========================================
cat("\n=== SAVING FILTERED OBJECT ===\n")
output_file <- file.path(output_dir, "merged_seurat_after_QC.rds")
saveRDS(merged_seurat_filtered, output_file)
cat("Saved filtered object to:", output_file, "\n")
cat("✓ Object ready for downstream analysis (normalization, scaling, etc.)\n")

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
cat("✓ Layers unified into single count matrix\n")
cat("✓ Metadata cleaned and organized\n")
cat("✓ Ribosomal genes (^RP[SL]) intentionally KEPT\n")

summary_table <- data.frame(
  Metric = c("Total Cells", "Cells Removed", "Cells Retained (%)", "Total Genes", 
             "Genes Removed", "Genes Retained (%)", "Median UMI/Cell", "Median Genes/Cell", 
             "Median Mito %", "Median Ribo %"),
  Before_Filtering = c(
    cells_before, "-", "100.00", genes_before, "-", "100.00",
    round(median(merged_seurat$nCount_RNA), 0),
    round(median(merged_seurat$nFeature_RNA), 0),
    round(median(merged_seurat$percent_mito), 2),
    round(median(merged_seurat$percent_ribo), 2)
  ),
  After_Filtering = c(
    cells_after, cells_removed, round(100 - removal_rate, 2),
    genes_after, genes_removed, round(100 - removal_rate_genes, 2),
    round(median(merged_seurat_filtered$nCount_RNA), 0),
    round(median(merged_seurat_filtered$nFeature_RNA), 0),
    round(median(merged_seurat_filtered$percent_mito), 2),
    round(median(merged_seurat_filtered$percent_ribo), 2)
  ),
  stringsAsFactors = FALSE
)

cat("\n=== FILTERING SUMMARY TABLE ===\n")
print(summary_table, row.names = FALSE)

gene_removal_details <- data.frame(
  Category = c("Mitochondrial (^MT-)", "Haemoglobin (^HB[^(P|E|S)])", "MALAT1",
               "Y-linked genes", "Low expression (< 5 cells)", "TOTAL (unique)"),
  Genes_Removed = c(mt_genes_removed, hb_genes_removed, malat1_removed, y_genes_removed,
                    length(low_expression_genes), length(all_genes_to_remove)),
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

excel_file <- file.path(output_dir, "QC_filtering_summary.xlsx")
write_xlsx(list("Summary" = summary_table, "Gene_Removal_Details" = gene_removal_details), excel_file)
cat("\n✓ Summary tables saved to:", excel_file, "\n")

cat("\nNEXT STEPS:\n")
cat("  → Normalization (SCTransform recommended for integrated data)\n")
cat("  → Scaling and variable gene selection\n")
cat("  → PCA\n")
cat("  → Batch correction if needed (Harmony, scVI)\n")
cat("  → Clustering and UMAP visualization\n")
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
