# Script SoupX - VERSIONE CORRETTA CON CLUSTERING DI CELL RANGER
library(SoupX)
library(Seurat)
library(Matrix)
library(dplyr)

# ============================================================================
# CONFIGURAZIONE
# ============================================================================

base_dir_1 <- "/home/fdann/Desktop/new_res"
base_dir_2 <- "/home/fdann/Desktop/wharf/fdann/fdann-sens2025518/pre-processing"

samples_dir1 <- c("P34304_1001", "P34304_1007", "P34304_1008", "P34304_1009")
samples_dir2 <- c("P33554_1001", "P33554_1002", "P33554_1003", "P33554_1004", 
                  "P34304_1002", "P34304_1003", "P34304_1004", "P34304_1005")

all_samples <- c(
  paste0(base_dir_1, "/", samples_dir1),
  paste0(base_dir_2, "/", samples_dir2)
)

sample_names <- c(samples_dir1, samples_dir2)

output_base_dir <- "/home/fdann/Desktop/proj/sc-res/0_soupx_analysis"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# FUNZIONE CORRETTA CON CLUSTERING DI CELL RANGER
# ============================================================================

apply_soupx_resolved <- function(sample_path, sample_name, output_dir) {
  tryCatch({
    cat("\n========================================\n")
    cat("Processando:", sample_name, "\n")
    
    # ========================================================================
    # VERIFICA STRUTTURA DIRECTORY
    # ========================================================================
    
    cat("📁 [1/8] Verifica struttura directory...\n")
    outs_path <- file.path(sample_path, "outs")
    if (!dir.exists(outs_path)) {
      stop("❌ Directory 'outs' non trovata in: ", sample_path)
    }
    
    filtered_path <- file.path(outs_path, "filtered_feature_bc_matrix")
    raw_path <- file.path(outs_path, "raw_feature_bc_matrix")
    analysis_path <- file.path(outs_path, "analysis")
    
    cat("   • outs_path:", outs_path, "\n")
    cat("   • filtered_path:", filtered_path, "\n")
    cat("   • raw_path:", raw_path, "\n")
    cat("   • analysis_path:", analysis_path, "\n")
    
    if (!dir.exists(filtered_path)) {
      stop("❌ Matrice FILTERED non trovata: ", filtered_path)
    }
    if (!dir.exists(analysis_path)) {
      stop("❌ Directory analysis non trovata: ", analysis_path)
    }
    cat("   ✅ Directory verificate\n")
    
    # ========================================================================
    # CARICAMENTO MATRICI
    # ========================================================================
    
    cat("📊 [2/8] Caricamento matrici...\n")
    
    # Carica matrice FILTERED (cellule reali)
    cat("   • Caricamento matrice FILTERED...\n")
    filtered_data <- Read10X(data.dir = filtered_path)
    cat("   ✅ Matrice FILTERED caricata\n")
    cat("     - Cellule:", ncol(filtered_data), "\n")
    cat("     - Geni:", nrow(filtered_data), "\n")
    
    # Carica matrice RAW (per soup background)
    cat("   • Caricamento matrice RAW...\n")
    raw_data <- Read10X(data.dir = raw_path)
    cat("   ✅ Matrice RAW caricata\n")
    cat("     - Barcodes totali:", ncol(raw_data), "\n")
    
    # ========================================================================
    # CARICAMENTO CLUSTERING DA CELL RANGER
    # ========================================================================
    
    cat("🔍 [3/8] Caricamento clustering...\n")
    
    # Prova a caricare il clustering graph-based (più comune)
    clustering_path <- file.path(analysis_path, "clustering", "gene_expression_graphclust", "clusters.csv")
    
    if (!file.exists(clustering_path)) {
      # Prova con kmeans come fallback
      clustering_path <- file.path(analysis_path, "clustering", "gene_expression_kmeans_10_clusters", "clusters.csv")
    }
    
    if (!file.exists(clustering_path)) {
      # Cerca qualsiasi file di clustering
      cluster_files <- list.files(file.path(analysis_path, "clustering"), 
                                 pattern = "clusters.csv", recursive = TRUE, full.names = TRUE)
      if (length(cluster_files) > 0) {
        clustering_path <- cluster_files[1]
        cat("   • Trovato file clustering alternativo:", clustering_path, "\n")
      } else {
        stop("❌ Nessun file di clustering trovato nella directory analysis")
      }
    }
    
    cat("   • Caricamento cluster da:", clustering_path, "\n")
    clusters_df <- read.csv(clustering_path)
    cat("   ✅ Cluster caricati\n")
    cat("     - Numero cluster:", length(unique(clusters_df$Cluster)), "\n")
    cat("     - Prime assegnazioni:", head(clusters_df$Cluster), "\n")
    
    # Prepara i cluster per SoupX
    clusters <- setNames(clusters_df$Cluster, clusters_df$Barcode)
    
    # Verifica che tutti i barcode nella matrice siano presenti nei cluster
    filtered_barcodes <- colnames(filtered_data)
    missing_barcodes <- setdiff(filtered_barcodes, names(clusters))
    
    if (length(missing_barcodes) > 0) {
      cat("   ⚠️  Barcodes mancanti nei cluster:", length(missing_barcodes), "\n")
      # Aggiungi cluster fittizio per i barcodes mancanti
      clusters[missing_barcodes] <- "Missing"
    }
    
    # Ordina i cluster per corrispondere all'ordine della matrice
    clusters <- clusters[filtered_barcodes]
    
    # ========================================================================
    # CREAZIONE SOUPCHANNEL E IMPOSTAZIONE CLUSTER
    # ========================================================================
    
    cat("🔧 [4/8] Creazione SoupChannel...\n")
    sc <- SoupChannel(tod = raw_data, toc = filtered_data)
    
    # Imposta i cluster
    cat("   • Impostazione cluster in SoupChannel...\n")
    sc <- setClusters(sc, clusters)
    cat("   ✅ SoupChannel creato con cluster\n")
    cat("     - Cellule nel channel:", ncol(sc$toc), "\n")
    cat("     - Cluster unici:", length(unique(sc$metaData$clusters)), "\n")
    
    # ========================================================================
    # ANALISI DATI PRIMA DELLA CORREZIONE
    # ========================================================================
    
    cat("📈 [5/8] Analisi dati pre-correzione...\n")
    cat("   • Statistiche matrice originale:\n")
    cat("     - Media conteggi per cellula:", round(mean(colSums(sc$toc)), 1), "\n")
    cat("     - Mediana conteggi per cellula:", round(median(colSums(sc$toc)), 1), "\n")
    
    # ========================================================================
    # STIMA CONTAMINAZIONE
    # ========================================================================
    
    cat("🎯 [6/8] Stima contaminazione...\n")
    cat("   • Esecuzione autoEstCont...\n")
    
    sc <- autoEstCont(sc, 
                      doPlot = FALSE,
                      forceAccept = TRUE)
    
    rho <- sc$fit$rhoEst
    cat("   ✅ Contaminazione stimata\n")
    cat("     - Rho globale:", round(rho, 4), "\n")
    cat("     - Percentuale:", round(rho * 100, 2), "%\n")
    
    # ========================================================================
    # APPLICAZIONE CORREZIONE
    # ========================================================================
    
    cat("🛠️  [7/8] Applicazione correzione SoupX...\n")
    cat("   • Esecuzione adjustCounts...\n")
    out <- adjustCounts(sc, 
                       roundToInt = TRUE,
                       method = "subtraction",
                       verbose = FALSE)
    
    cat("   ✅ Correzione applicata\n")
    cat("     - Cellule dopo correzione:", ncol(out), "\n")
    
    # Verifica coerenza dimensionale
    if (ncol(out) != ncol(filtered_data)) {
      cat("   ❌ ALLERTA: Numero cellule cambiato!\n")
    } else {
      cat("   ✅ Numero cellule coerente\n")
    }
    
    # ========================================================================
    # SALVATAGGIO DATI
    # ========================================================================
    
    cat("💾 [8/8] Salvataggio dati...\n")
    output_path <- file.path(output_dir, paste0(sample_name, "_soupX_corrected"))
    
    # Pulizia directory esistente
    if (dir.exists(output_path)) {
      unlink(output_path, recursive = TRUE, force = TRUE)
      Sys.sleep(1)
    }
    
    # Creazione directory
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    
    # Salvataggio manuale
    cat("   • Salvataggio file...\n")
    writeMM(out, file = file.path(output_path, "matrix.mtx"))
    
    write.table(
      colnames(out),
      file = file.path(output_path, "barcodes.tsv"),
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    features_df <- data.frame(
      gene_id = rownames(out),
      gene_name = rownames(out),
      feature_type = "Gene Expression"
    )
    write.table(
      features_df,
      file = file.path(output_path, "features.tsv"),
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    # Compressione
    cat("   • Compressione file...\n")
    system(paste("gzip -f", file.path(output_path, "matrix.mtx")))
    system(paste("gzip -f", file.path(output_path, "barcodes.tsv")))
    system(paste("gzip -f", file.path(output_path, "features.tsv")))
    
    # Verifica finale
    files <- list.files(output_path, pattern = "\\.gz$")
    cat("   ✅ File creati:", paste(files, collapse = ", "), "\n")
    cat("✅", sample_name, "COMPLETATO\n")
    
    return(list(
      success = TRUE,
      sample_name = sample_name,
      output_path = output_path,
      contamination = rho,
      cells = ncol(out)
    ))
    
  }, error = function(e) {
    cat("❌ ERRORE in", sample_name, ":\n")
    cat("   ", e$message, "\n")
    return(list(
      success = FALSE,
      sample_name = sample_name,
      error = e$message
    ))
  })
}

# ============================================================================
# ESECUZIONE PRINCIPALE
# ============================================================================

cat("INIZIO PROCESSAMENTO SOUPX - CON CLUSTERING CELL RANGER\n")
cat("Numero campioni:", length(all_samples), "\n")
cat("========================================\n")

results <- list()

for (i in seq_along(all_samples)) {
  sample_path <- all_samples[i]
  sample_name <- sample_names[i]
  
  cat("\n", i, "/", length(all_samples), "-", sample_name, "\n")
  
  result <- apply_soupx_resolved(sample_path, sample_name, output_base_dir)
  results[[sample_name]] <- result
}

# ============================================================================
# RESOCONTO E CREAZIONE SEURAT
# ============================================================================

successful_samples <- names(results)[sapply(results, function(x) x$success)]

cat("\n\n========================================\n")
cat("RESOCONTO FINALE\n")
cat("========================================\n")
cat("Campioni processati con successo:", length(successful_samples), "/", length(all_samples), "\n\n")

if (length(successful_samples) > 0) {
  cat("DETTAGLI:\n")
  for (name in successful_samples) {
    result <- results[[name]]
    cat(sprintf("✓ %s: %.2f%% contaminazione, %d cellule\n", 
                name, result$contamination * 100, result$cells))
  }
  
  # CREAZIONE SEURAT SENZA FILTRI
  cat("\nCreazione oggetto Seurat merged...\n")
  seurat_list <- list()
  
  for (sample_name in successful_samples) {
    output_path <- results[[sample_name]]$output_path
    
    # Verifica che la directory contenga file
    files <- list.files(output_path)
    if (length(files) >= 3) {
      data <- Read10X(data.dir = output_path)
      
      # CREAZIONE SEURAT SENZA FILTRI
      seurat_obj <- CreateSeuratObject(
        counts = data,
        project = sample_name
      )
      
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      seurat_obj[["sample"]] <- sample_name
      seurat_obj[["soupx_contamination"]] <- results[[sample_name]]$contamination
      
      seurat_list[[sample_name]] <- seurat_obj
      cat("✓", sample_name, "-", ncol(seurat_obj), "cellule\n")
    } else {
      cat("✗", sample_name, "- directory incompleta\n")
    }
  }
  
  # MERGE
  if (length(seurat_list) > 0) {
    if (length(seurat_list) > 1) {
      merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], 
                            add.cell.ids = names(seurat_list))
    } else {
      merged_seurat <- seurat_list[[1]]
    }
    
    save_path <- file.path(output_base_dir, "merged_seurat_soupX_corrected.rds")
    saveRDS(merged_seurat, save_path)
    cat("\n✅ OGGETTO SEURAT SALVATO:", save_path, "\n")
    cat("   Cellule totali:", ncol(merged_seurat), "\n")
    cat("   Geni totali:", nrow(merged_seurat), "\n")
  }
} else {
  cat("\n❌ Nessun campione processato con successo\n")
}
