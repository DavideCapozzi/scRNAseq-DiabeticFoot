# ==============================================================================
# SCRIPT 00: DOWNLOAD E VERIFICA CAMPIONI GSE165816
# Dataset: Kim et al., Nature 2023 - Diabetic Foot Ulcer scRNA-seq
#
# GRUPPI:
#   - Diabetic Foot Ulcer Healing    (9 pazienti, 11 campioni)
#   - Diabetic Foot Ulcer Non-Healing (5 pazienti,  6 campioni)
#   - Diabetes Control               (6 pazienti,  8 campioni)  ← NUOVO
#
# CORREZIONE RISPETTO ALLA VERSIONE PRECEDENTE:
#   - Aggiunto gruppo Diabetes (prima assente)
#   - Aggiunta verifica post-download per ogni file
#   - Gestione robusta degli ID con suffisso "A" (es. G1A, G3A, G2A)
#   - Report di verifica finale
# ==============================================================================

rm(list = ls())

library(dplyr)
library(stringr)
library(GEOquery)

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/nature_paper_validation"
data_dir <- file.path(base_dir, "scRNA_data_new")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. DEFINIZIONE CAMPIONI (hardcoded dalla tabella del paper + foglio Excel)
#    NOTA: Gli ID con "&" indicano campioni multiplexati dallo stesso paziente
#          (es. "G4 & G23" = stesso paziente, due file separati su GEO)
#          Gli ID con "A" (G1A, G2A, G3A) sono skin biopsies di controllo
#          adiacenti nei pazienti diabetici
# ==============================================================================
sample_groups <- list(

  "Healing" = c(
    "G4", "G23",   # stesso paziente
    "G2",
    "G15",
    "G7", "G8",    # stesso paziente
    "G49",
    "G42",
    "G45"
  ),

  "Non-Healing" = c(
    "G6",
    "G9",
    "G33", "G34",  # stesso paziente
    "G39"
  ),

  "Diabetes" = c(    # ← GRUPPO CONTROLLO - PRIMA MANCANTE
    "G2A",
    "G46",
    "G38",
    "G41",
    "G1A", "G3A",  # stesso paziente
    "G3", "G5"     # stesso paziente
  )
)

# Crea un data frame flat per il mapping
metadata_samples <- do.call(rbind, lapply(names(sample_groups), function(grp) {
  data.frame(
    Sample_ID = sample_groups[[grp]],
    Condition = grp,
    stringsAsFactors = FALSE
  )
}))

cat("=== CAMPIONI DA SCARICARE ===\n")
cat(sprintf("  Healing    : %d campioni\n", sum(metadata_samples$Condition == "Healing")))
cat(sprintf("  Non-Healing: %d campioni\n", sum(metadata_samples$Condition == "Non-Healing")))
cat(sprintf("  Diabetes   : %d campioni\n", sum(metadata_samples$Condition == "Diabetes")))
cat(sprintf("  TOTALE     : %d campioni\n", nrow(metadata_samples)))
print(metadata_samples)

# Salva tabella campioni per uso negli script successivi
write.csv(metadata_samples,
          file.path(base_dir, "sample_metadata.csv"),
          row.names = FALSE)
cat("\nTabella metadati salvata in: sample_metadata.csv\n")

# ==============================================================================
# 2. CONNESSIONE A GEO E MAPPING GSM ↔ SAMPLE_ID
# ==============================================================================
cat("\n--- Connessione a GEO (GSE165816) ---\n")
cat("Questo può richiedere 1-2 minuti...\n")

gse_metadata <- tryCatch(
  getGEO("GSE165816", GSEMatrix = FALSE),
  error = function(e) {
    stop("Impossibile connettersi a GEO. Verifica la connessione internet.\nErrore: ", e$message)
  }
)

# Estrai titoli di tutti i GSM
gsm_titles <- unlist(lapply(GSMList(gse_metadata), function(x) Meta(x)$title))
cat(sprintf("GSM totali nel dataset: %d\n", length(gsm_titles)))

# Pattern matching robusto per tutti i Sample ID
# NOTA: per G1A, G3A, G2A usiamo \\b per evitare match parziali
#       Alcuni ID potrebbero apparire come "G1 adjacent" o "G1A" nel titolo GEO
all_ids <- metadata_samples$Sample_ID


# Costruisci pattern che gestisce sia "G2A" che possibili varianti
build_pattern <- function(id) {
  # Trasforma "G2A" in pattern che cattura anche varianti
  if (grepl("[A-Z]$", id) && nchar(id) > 2) {
    # ID con suffisso lettera (es G1A): cerca sia "G1A" che "G1 A" che "G1-A"
    base <- sub("[A-Z]$", "", id)
    suffix <- sub(".*([A-Z])$", "\\1", id)
    return(paste0("\\b", id, "\\b|", base, "\\s*", suffix, "\\b"))
  } else {
    return(paste0("\\b", id, "\\b"))
  }
}

patterns_per_id <- sapply(all_ids, build_pattern)

# Filtra GSM: deve matchare un ID E essere "Skin/Foot" tissue
gsm_mapping <- data.frame(
  GSM_ID      = character(),
  GSM_Title   = character(),
  Sample_ID   = character(),
  Condition   = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(all_ids)) {
  id  <- all_ids[i]
  pat <- patterns_per_id[[id]]
  
  # Match su titolo (case-insensitive per sicurezza)
  matched <- names(gsm_titles)[
    grepl(pat, gsm_titles, ignore.case = FALSE) &
    grepl("Skin|Foot|skin|foot", gsm_titles, ignore.case = TRUE)
  ]
  
  if (length(matched) == 0) {
    cat(sprintf("  [WARN] Nessun GSM trovato per Sample ID: %s\n", id))
  } else {
    for (gsm in matched) {
      condition <- metadata_samples$Condition[metadata_samples$Sample_ID == id]
      gsm_mapping <- rbind(gsm_mapping, data.frame(
        GSM_ID    = gsm,
        GSM_Title = gsm_titles[[gsm]],
        Sample_ID = id,
        Condition = condition,
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  [OK] %s → %s (%s)\n", id, gsm, gsm_titles[[gsm]]))
    }
  }
}

# Rimuovi duplicati (un GSM mappato a più ID → prendi il primo)
gsm_mapping <- gsm_mapping[!duplicated(gsm_mapping$GSM_ID), ]

cat(sprintf("\nGSM trovati: %d / %d richiesti\n",
            nrow(gsm_mapping), nrow(metadata_samples)))

# Salva mapping
write.csv(gsm_mapping, file.path(base_dir, "gsm_mapping.csv"), row.names = FALSE)

# ==============================================================================
# 3. DOWNLOAD DEI FILE SUPPLEMENTARI
# ==============================================================================
cat("\n--- Avvio download ---\n")
cat(sprintf("Download di %d campioni in: %s\n\n", nrow(gsm_mapping), data_dir))

download_log <- data.frame(
  GSM_ID    = character(),
  Sample_ID = character(),
  Condition = character(),
  Status    = character(),
  Files     = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(gsm_mapping))) {
  gsm_id    <- gsm_mapping$GSM_ID[i]
  sample_id <- gsm_mapping$Sample_ID[i]
  condition <- gsm_mapping$Condition[i]
  
  cat(sprintf("[%d/%d] Scaricando %s (%s - %s)...\n",
              i, nrow(gsm_mapping), gsm_id, sample_id, condition))
  
  gsm_dir <- file.path(data_dir, gsm_id)
  
  status <- tryCatch({
    getGEOSuppFiles(gsm_id, makeDirectory = TRUE, baseDir = data_dir)
    "SUCCESS"
  }, error = function(e) {
    cat(sprintf("  [ERROR] %s\n", e$message))
    paste0("ERROR: ", e$message)
  })
  
  # Registra file scaricati
  files_found <- if (dir.exists(gsm_dir)) {
    paste(list.files(gsm_dir), collapse = "; ")
  } else {
    "NO DIR"
  }
  
  download_log <- rbind(download_log, data.frame(
    GSM_ID    = gsm_id,
    Sample_ID = sample_id,
    Condition = condition,
    Status    = status,
    Files     = files_found,
    stringsAsFactors = FALSE
  ))
  
  Sys.sleep(0.5)  # pausa breve per non sovraccaricare GEO
}

# ==============================================================================
# 4. VERIFICA POST-DOWNLOAD
# ==============================================================================
cat("\n=== VERIFICA POST-DOWNLOAD ===\n")

verify_download <- function(gsm_id, data_dir) {
  gsm_dir <- file.path(data_dir, gsm_id)
  
  if (!dir.exists(gsm_dir)) return(list(ok = FALSE, reason = "Directory non trovata"))
  
  files <- list.files(gsm_dir, full.names = TRUE,
                       pattern = "\\.(gz|csv|tsv|txt|mtx)$",
                       recursive = TRUE)
  
  if (length(files) == 0) return(list(ok = FALSE, reason = "Nessun file trovato"))
  
  # Verifica dimensione > 1KB
  sizes <- file.info(files)$size
  if (any(sizes < 1024)) {
    small <- files[sizes < 1024]
    return(list(ok = FALSE,
                reason = paste("File troppo piccolo (<1KB):", basename(small[1]))))
  }
  
  # Verifica che almeno un .gz sia leggibile
  gz_files <- files[grepl("\\.gz$", files)]
  if (length(gz_files) > 0) {
    readable <- tryCatch({
      con <- gzfile(gz_files[1], "r")
      first_line <- readLines(con, n = 1)
      close(con)
      nchar(first_line) > 0
    }, error = function(e) FALSE)
    
    if (!readable) return(list(ok = FALSE, reason = "File .gz corrotto"))
  }
  
  return(list(ok = TRUE, reason = paste(length(files), "file OK")))
}

cat("\nVerifica integrità per ogni campione:\n")
verification_results <- lapply(gsm_mapping$GSM_ID, function(gsm) {
  result <- verify_download(gsm, data_dir)
  status_icon <- if (result$ok) "[✓]" else "[✗]"
  idx <- which(gsm_mapping$GSM_ID == gsm)
  cat(sprintf("  %s %s (%s) - %s\n",
              status_icon, gsm,
              gsm_mapping$Sample_ID[idx],
              result$reason))
  data.frame(
    GSM_ID    = gsm,
    Sample_ID = gsm_mapping$Sample_ID[idx],
    Condition = gsm_mapping$Condition[idx],
    Verified  = result$ok,
    Notes     = result$reason,
    stringsAsFactors = FALSE
  )
})

verification_df <- do.call(rbind, verification_results)

n_ok  <- sum(verification_df$Verified)
n_bad <- sum(!verification_df$Verified)

cat(sprintf("\n=== SOMMARIO VERIFICA ===\n"))
cat(sprintf("  Download OK  : %d / %d\n", n_ok, nrow(verification_df)))
cat(sprintf("  Download FAIL: %d / %d\n", n_bad, nrow(verification_df)))

if (n_bad > 0) {
  cat("\n[ATTENZIONE] Campioni con problemi:\n")
  bad_samples <- verification_df[!verification_df$Verified, ]
  print(bad_samples[, c("GSM_ID", "Sample_ID", "Condition", "Notes")])
  cat("\nRiprovare a scaricare questi campioni manualmente o verificare\n")
  cat("il mapping GSM in:", file.path(base_dir, "gsm_mapping.csv"), "\n")
}

# Salva report verifica
write.csv(verification_df,
          file.path(base_dir, "download_verification_report.csv"),
          row.names = FALSE)

cat("\n=== DOWNLOAD COMPLETATO ===\n")
cat("Mapping GSM salvato      :", file.path(base_dir, "gsm_mapping.csv"), "\n")
cat("Report verifica salvato  :", file.path(base_dir, "download_verification_report.csv"), "\n")
cat("Metadati campioni salvati:", file.path(base_dir, "sample_metadata.csv"), "\n")
cat("\nProcedi con: 01_qc_preprocessing.R\n")

