rm(list = ls())

library(Seurat)
library(ggplot2)
library(biomaRt)
library(xlsx)

# ============================================================================
# CONFIGURAZIONE INIZIALE
# ============================================================================

# Percorsi input/output
#input_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/1_pre_filtering/merged_seurat_before_QC.rds"
# output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/1_pre_filtering/qc_violin_plots_samples"
# output_seurat <- "G:/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/1_pre_filtering/merged_seurat_withAllStats.rds"

#input_file <- "~/Desktop/proj/sc-res/2_new_analysis/1_pre_processing/merged_seurat_after_doublet_QC.rds"
# output_dir <- "~/Desktop/proj/sc-res/2_new_analysis/1_pre_processing/qc_violin_plots_samples"
# output_seurat <- "~/Desktop/proj/sc-res/2_new_analysis/1_pre_processing/merged_seurat_withAllStats.rds"

input_file <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/2_pre_filtering/merged_seurat_after_doublet_QC.rds"
output_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/2_pre_filtering/qc_violin_plots_samples"
patient_info_file <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/00_new_analysis/info_patients/patient_info.xlsx"
output_seurat <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/2_pre_filtering/merged_seurat_withAllStats.rds"
  
# Colonna per raggruppamento nei plot
sample_column <- "orig.ident"

# Colonna per colorazione nei plot
color_column <- "condition"  

# Colori personalizzati per condizione
condition_colors <- c(
  "not_healed" = "#fdae61",
  "healed" = "#7fbc41"
)

sample_order <- c(
  "P33554_1003", "P33554_1004", "P34304_1002", "P34304_1003", 
  "P34304_1004", "P34304_1008", "P34304_1009", "P33554_1001", 
  "P33554_1002", "P34304_1001", "P34304_1005", "P34304_1007"
)

# Assay da analizzare
assays_to_analyze <- c("RNA")

calculate_stats <- TRUE  # Cambia in FALSE per usare stats già calcolate

# SOGLIE QC (facilmente modificabili)
thresholds <- list(
  nFeature_RNA = list(min = 200, max = NULL),  # min: 200, max: NULL (no upper limit)
  nCount_RNA = list(min = 1000, max = NULL),    # nessuna soglia di default
  percent_mito = list(min = NULL, max = 20),    # max: 20%
  percent_ribo = list(min = 5, max = NULL),     # min: 5%
  percent_hb = list(min = NULL, max = NULL),
  percent_malat1 = list(min = NULL, max = NULL),
  percent_Y_genes = list(min = NULL, max = NULL),
  percent_XIST = list(min = NULL, max = NULL)
)

# ============================================================================
# STEP 1: CARICAMENTO DATI E CALCOLO STATISTICHE
# ============================================================================

message("=== STEP 1: CARICAMENTO E CALCOLO STATISTICHE ===")
seurat_obj <- readRDS(input_file)

# Carica la tabella patient_info
message("Caricamento patient_info.xlsx...")
patient_info <- read.xlsx(patient_info_file, sheetIndex  = "Data")

# Verifica la struttura della tabella
message("Struttura patient_info:")
print(head(patient_info))
message(paste("Campioni in patient_info:", nrow(patient_info)))

# Verifica i campioni presenti nell'oggetto Seurat
seurat_samples <- unique(seurat_obj@meta.data[[sample_column]])
message(paste("Campioni in Seurat object:", length(seurat_samples)))
message("Campioni Seurat:", paste(seurat_samples, collapse = ", "))

# Crea una mappa Fastq_File_ID -> outcome
sample_to_condition <- setNames(patient_info$outcome, patient_info$Fastq_File_ID)
message("Mappa campione -> condition:")
print(sample_to_condition)

# Verifica che tutti i campioni in Seurat siano presenti nella patient_info
missing_samples <- setdiff(seurat_samples, names(sample_to_condition))
if (length(missing_samples) > 0) {
  warning("ATTENZIONE: Campioni in Seurat non trovati in patient_info: ", 
          paste(missing_samples, collapse = ", "))
} else {
  message("✓ Tutti i campioni Seurat trovati in patient_info")
}

# Aggiungi la colonna condition ai metadati
message("Aggiunta colonna 'condition' ai metadati...")
seurat_obj@meta.data$condition <- sample_to_condition[as.character(seurat_obj@meta.data[[sample_column]])]

# Verifica l'aggiunta
message("Verifica metadati aggiunti:")
print(table(seurat_obj@meta.data[[sample_column]], seurat_obj@meta.data$condition, useNA = "ifany"))

# ----------------------------------------------------------------------------
# 1.1 Funzione per calcolare statistiche QC standard
# ----------------------------------------------------------------------------
calculate_qc_stats <- function(seurat_obj, assay_name) {
  
  DefaultAssay(seurat_obj) <- assay_name
  message(paste("Processando assay:", assay_name))
  
  current_genes <- rownames(seurat_obj)
  stats_calculated <- c("nFeature_RNA", "nCount_RNA")
  
  # Geni mitocondriali
  mt_genes <- grep("^MT-", current_genes, value = TRUE, ignore.case = TRUE)
  if (length(mt_genes) > 0) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", 
                                       col.name = "percent_mito")
    stats_calculated <- c(stats_calculated, "percent_mito")
    message(paste("  ✓ Geni mitocondriali:", length(mt_genes)))
  } else {
    message("  ⚠ Nessun gene mitocondriale trovato")
  }
  
  # Geni ribosomiali
  rp_genes <- grep("^RP[SL]", current_genes, value = TRUE, ignore.case = TRUE)
  if (length(rp_genes) > 0) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]", 
                                       col.name = "percent_ribo")
    stats_calculated <- c(stats_calculated, "percent_ribo")
    message(paste("  ✓ Geni ribosomiali:", length(rp_genes)))
  } else {
    message("  ⚠ Nessun gene ribosomiale trovato")
  }
  
  # Geni emoglobina
  hb_genes <- grep("^HB[^(P|E|S)]", current_genes, value = TRUE, ignore.case = TRUE)
  if (length(hb_genes) > 0) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P|E|S)]", 
                                       col.name = "percent_hb")
    stats_calculated <- c(stats_calculated, "percent_hb")
    message(paste("  ✓ Geni emoglobina:", length(hb_genes)))
  } else {
    message("  ⚠ Nessun gene emoglobina trovato")
  }
  
  # MALAT1
  if ("MALAT1" %in% current_genes) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MALAT1$", 
                                       col.name = "percent_malat1")
    stats_calculated <- c(stats_calculated, "percent_malat1")
    message("  ✓ MALAT1 trovato")
  } else {
    message("  ⚠ MALAT1 non trovato")
  }
  
  return(list(seurat_obj = seurat_obj, stats = stats_calculated))
}

# ----------------------------------------------------------------------------
# 1.2 Funzione per calcolare statistiche geni sessuali (solo per assay RNA)
# ----------------------------------------------------------------------------
calculate_sex_genes_stats <- function(seurat_obj) {
  
  message("Calcolando statistiche geni sessuali...")
  current_genes <- rownames(seurat_obj)
  stats_calculated <- c()
  
  # Connessione a Ensembl
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Geni cromosoma Y
  genes_Y <- getBM(attributes = c("hgnc_symbol", "chromosome_name"),
                   filters = "chromosome_name",
                   values = "Y",
                   mart = mart)
  
  genes_X <- getBM(attributes = c("hgnc_symbol", "chromosome_name"),
                   filters = "chromosome_name",
                   values = "X",
                   mart = mart)
  
  genes_Y_only <- setdiff(genes_Y$hgnc_symbol, genes_X$hgnc_symbol)
  genes_Y_in_data <- genes_Y_only[genes_Y_only %in% current_genes]
  
  if (length(genes_Y_in_data) > 0) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, features = genes_Y_in_data, 
                                       col.name = "percent_Y_genes")
    stats_calculated <- c(stats_calculated, "percent_Y_genes")
    message(paste("  ✓ Geni Y trovati:", length(genes_Y_in_data)))
  } else {
    message("  ⚠ Nessun gene Y trovato")
  }
  
  # Gene XIST
  if ("XIST" %in% current_genes) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, features = "XIST", 
                                       col.name = "percent_XIST")
    stats_calculated <- c(stats_calculated, "percent_XIST")
    message("  ✓ Gene XIST trovato")
  } else {
    message("  ⚠ Gene XIST non trovato")
  }
  
  return(list(seurat_obj = seurat_obj, stats = stats_calculated))
}

# ----------------------------------------------------------------------------
# 1.3 Calcolo statistiche per tutti gli assay
# ----------------------------------------------------------------------------
all_features <- list()

if (calculate_stats) {
  # CALCOLA LE STATISTICHE
  message("\n→ Modalità: CALCOLO NUOVE STATISTICHE")
  
  for (assay in assays_to_analyze) {
    if (assay %in% names(seurat_obj@assays)) {
      result <- calculate_qc_stats(seurat_obj, assay)
      seurat_obj <- result$seurat_obj
      all_features[[assay]] <- result$stats
    } else {
      message(paste("  ⚠ Assay", assay, "non trovato nell'oggetto Seurat"))
    }
  }
  
  # Statistiche geni sessuali (solo per RNA)
  if ("RNA" %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- "RNA"
    sex_result <- calculate_sex_genes_stats(seurat_obj)
    seurat_obj <- sex_result$seurat_obj
    all_features[["RNA"]] <- c(all_features[["RNA"]], sex_result$stats)
  }
  
} else {
  # USA STATISTICHE GIÀ PRESENTI
  message("\n→ Modalità: USO STATISTICHE ESISTENTI")
  
  # Lista di tutte le possibili stats QC
  possible_stats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", 
                      "percent_ribo", "percent_hb", "percent_malat1",
                      "percent_Y_genes", "percent_XIST")
  
  # Rileva quali stats sono presenti nei metadata
  existing_stats <- possible_stats[possible_stats %in% colnames(seurat_obj@meta.data)]
  
  if (length(existing_stats) == 0) {
    stop("ERRORE: Nessuna statistica QC trovata nei metadata. Imposta calculate_stats = TRUE")
  }
  
  # Assegna le stats trovate agli assay specificati
  for (assay in assays_to_analyze) {
    if (assay %in% names(seurat_obj@assays)) {
      all_features[[assay]] <- existing_stats
      message(paste("  ✓ Assay", assay, "- stats trovate:", 
                    paste(existing_stats, collapse = ", ")))
    }
  }
}

# ============================================================================
# STEP 2: SALVATAGGIO OGGETTO SEURAT
# ============================================================================

message("\n=== STEP 2: SALVATAGGIO OGGETTO SEURAT ===")
if (calculate_stats) {
  saveRDS(seurat_obj, output_seurat)
  message(paste("✓ Oggetto Seurat salvato:", output_seurat))
} else {
  message("⊘ Salvataggio saltato (stats non ricalcolate)")
}

# ============================================================================
# STEP 3: CREAZIONE VIOLIN PLOTS
# ============================================================================

message("\n=== STEP 3: CREAZIONE VIOLIN PLOTS ===")

# Crea directory output se non esiste
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("✓ Directory creata:", output_dir))
}

# Verifica colonna campioni
if (!sample_column %in% colnames(seurat_obj@meta.data)) {
  stop(paste("Errore: La colonna", sample_column, "non esiste nei metadata!"))
}

seurat_obj@meta.data[[sample_column]] <- factor(
  seurat_obj@meta.data[[sample_column]], 
  levels = sample_order
)

# ----------------------------------------------------------------------------
# 3.1 Funzione per creare violin plot con threshold
# ----------------------------------------------------------------------------
create_violin_with_threshold <- function(seurat_obj, feature, assay_name, 
                                         sample_col, thresh_list, color_col, cond_colors) {
  
  # Crea plot base
  p <- VlnPlot(seurat_obj, 
               features = feature, 
               group.by = sample_col,
               split.by = color_column,
               pt.size = 0) +
    scale_fill_manual(values = condition_colors) +
    ggtitle(paste0(feature, " - ", assay_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.x = element_blank()
    )
  
  # Aggiungi linee di threshold se presenti
  if (feature %in% names(thresh_list)) {
    thresh <- thresh_list[[feature]]
    
    if (!is.null(thresh$min)) {
      p <- p + geom_hline(yintercept = thresh$min, 
                          linetype = "dashed", 
                          color = "red", 
                          linewidth = 0.8)
    }
    
    if (!is.null(thresh$max)) {
      p <- p + geom_hline(yintercept = thresh$max, 
                          linetype = "dashed", 
                          color = "red", 
                          linewidth = 0.8)
    }
  }
  
  return(p)
}

# ----------------------------------------------------------------------------
# 3.2 Creazione e salvataggio plot per ogni feature
# ----------------------------------------------------------------------------
for (assay in names(all_features)) {
  
  message(paste("\nCreando plot per assay:", assay))
  DefaultAssay(seurat_obj) <- assay
  
  for (feature in all_features[[assay]]) {
    
    tryCatch({
      
      # Crea plot
      p <- create_violin_with_threshold(seurat_obj, feature, assay, 
                                        sample_column, thresholds, 
                                        color_column, condition_colors)
      
      # Salva come PDF (useCairoPDF per file più leggeri)
      file_name <- file.path(output_dir, 
                             paste0(assay, "_", feature, "_by_sample.pdf"))
      
      ggsave(file_name, p, 
             width = 10, height = 6, 
             device = "pdf",
             dpi = 150)  # DPI ridotto per file più leggeri
      
      message(paste("  ✓ Salvato:", basename(file_name)))
      
    }, error = function(e) {
      message(paste("  ✗ Errore nel creare plot per", feature, ":", e$message))
    })
  }
}

# ============================================================================
# RIEPILOGO FINALE
# ============================================================================

message("\n=== ANALISI COMPLETATA ===")
message(paste("✓ Oggetto Seurat con tutte le statistiche:", output_seurat))
message(paste("✓ Violin plots salvati in:", output_dir))
message("\nStatistiche calcolate per assay:")
for (assay in names(all_features)) {
  message(paste("  -", assay, ":", paste(all_features[[assay]], collapse = ", ")))
}
message("\nSoglie applicate nei plot:")
for (feat in names(thresholds)) {
  thresh <- thresholds[[feat]]
  if (!is.null(thresh$min) || !is.null(thresh$max)) {
    min_val <- if (is.null(thresh$min)) "none" else thresh$min
    max_val <- if (is.null(thresh$max)) "none" else thresh$max
    message(paste("  -", feat, ": min =", min_val, ", max =", max_val))
  }
}

