# ==============================================================================
# SCRIPT: Integrazione Proteomica vs Single Cell (Final + BioMart Integration)
# DESCRIZIONE: 
# Identifica i geni tramite:
# 1. Match Esatto
# 2. Alias Manuali
# 3. HGNC Helper (Simboli ufficiali/storici)
# 4. Org.Hs.eg.db (Entrez ID)
# 5. BioMart (Ensembl Synonyms - NUOVO LIVELLO)
# ==============================================================================

# 1. LIBRERIE ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(HGNChelper)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(biomaRt) # Aggiunta libreria BioMart
})

# 2. PATH E FILE DI INPUT ------------------------------------------------------
root_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE"

# File di input (Risultati DE Single Cell)
deg_file_path <- file.path(root_path, "01_second_new_analysis/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusters.xlsx")

# File di output (Aggiornato nome file)
output_file_path <- file.path(root_path, "01_second_new_analysis/6_differential_expression_analysis/results/DEGs_Genes_of_Interest_BioMart_FINAL.xlsx")

# 3. LISTE GENI (HARDCODED) ----------------------------------------------------
visit1_genes <- c("IL-15RA", "IL-12B", "TNF", "CCL3", "CCL4", "MCP-4", "IL-18R1")
visit2_genes <- c("SPARCL1", "FCN2", "VASN", "CHL1", "CCL18", "SAA4", "F11", "SELL", "F7", "LTBP2", "CNDP1", "CRTAC1", "ICAM3")

# Unione liste per query BioMart
all_input_genes <- unique(c(visit1_genes, visit2_genes))

# 4. DIZIONARIO MANUALE (Alias specifici) --------------------------------------
manual_alias_map <- list(
  "IL-15RA" = c("IL15RA", "CD215"),
  "IL-12B"  = c("IL12B", "NFSF2", "CLMF"),
  "MCP-4"   = c("CCL13", "NCC-1", "MCP4"),
  "IL-18R1" = c("IL18R1", "CD218A", "IL1RRP"),
  "F7"      = c("FVII", "SPCA"),
  "F11"     = c("FXI", "PTA"),
  "SELL"    = c("CD62L", "L-Selectin"),
  "ICAM3"   = c("CD50")
)

# 5. PRE-FETCHING BIOMART (NUOVO BLOCCO) ---------------------------------------
# Scarichiamo i sinonimi PRIMA di iniziare il loop per velocità e stabilità
cat("\nConfigurazione BioMart...\n")
biomart_lookup <- list() # Lista vuota di default

tryCatch({
  # Setup connessione
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Scarica dati solo per i nostri geni
  cat("Scaricamento sinonimi da BioMart per", length(all_input_genes), "geni...\n")
  bm_data <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
                   filters = "hgnc_symbol",
                   values = all_input_genes,
                   mart = mart)
  
  # Rimuovi sinonimi vuoti e organizza in lista
  bm_data <- bm_data %>% filter(external_synonym != "")
  
  for(g in all_input_genes) {
    syns <- bm_data %>% filter(hgnc_symbol == g) %>% pull(external_synonym)
    if(length(syns) > 0) {
      biomart_lookup[[g]] <- unique(syns)
    }
  }
  cat("BioMart: Sinonimi recuperati con successo.\n")
  
}, error = function(e) {
  warning("ATTENZIONE: Impossibile connettersi a BioMart. Lo script continuerà senza questo passaggio.")
  cat("Dettaglio errore BioMart:", e$message, "\n")
})

# 6. CARICAMENTO E PREPARAZIONE DATI (Safe Mode) -------------------------------
cat("Caricamento file DE Single Cell...\n")

if(!file.exists(deg_file_path)) stop("Errore critico: Il file di input non esiste.")

# Lettura file
sc_de_data <- read_excel(deg_file_path, sheet = 1)

# Standardizzazione nomi colonne
colnames(sc_de_data) <- tolower(colnames(sc_de_data))

# Identificazione Colonna Gene
if ("gene" %in% colnames(sc_de_data)) {
  cat("Colonna 'gene' trovata correttamente.\n")
} else if ("...1" %in% colnames(sc_de_data)) {
  cat("Uso la prima colonna ('...1') come gene.\n")
  sc_de_data <- sc_de_data %>% dplyr::rename(gene = `...1`)
} else {
  stop("Errore: Impossibile trovare una colonna valida per i geni.")
}

# Verifica Colonna Cluster
if (!"cluster" %in% colnames(sc_de_data)) {
  stop("Errore: La colonna 'cluster' non è stata trovata nel file Excel.")
}

# Creazione universo geni
sc_genes_universe <- unique(sc_de_data$gene)
cat("Dataset pronto. Analisi su", length(sc_genes_universe), "geni unici.\n")


# 7. FUNZIONI DI MATCHING (Aggiornata con BioMart) -----------------------------

clean_symbol <- function(s) {
  toupper(trimws(s))
}

find_gene_match <- function(input_gene, sc_universe) {
  clean_input <- clean_symbol(input_gene)
  
  # 1. Match Diretto
  if (clean_input %in% sc_universe) return(list(matched_gene = clean_input, type = "Direct"))
  
  # 2. Match Manuale
  if (input_gene %in% names(manual_alias_map)) {
    candidates <- manual_alias_map[[input_gene]]
    match <- intersect(candidates, sc_universe)
    if (length(match) > 0) return(list(matched_gene = match[1], type = "Manual_Alias"))
  }
  
  # 3. HGNC Helper (Ufficiale)
  hgnc_check <- suppressWarnings(checkGeneSymbols(clean_input, unmapped.as.na = FALSE))
  if (hgnc_check$Approved == FALSE && !is.na(hgnc_check$Suggested.Symbol)) {
    suggested <- strsplit(hgnc_check$Suggested.Symbol, " /// ")[[1]]
    match <- intersect(suggested, sc_universe)
    if (length(match) > 0) return(list(matched_gene = match[1], type = "HGNC_Correction"))
  }
  
  # 4. Org.Hs.eg.db (Database locale Entrez)
  tryCatch({
    query <- if(hgnc_check$Approved == FALSE && !is.na(hgnc_check$Suggested.Symbol)) hgnc_check$Suggested.Symbol[1] else clean_input
    entrez <- suppressMessages(mapIds(org.Hs.eg.db, keys = query, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
    if (!is.na(entrez)) {
      aliases <- suppressMessages(select(org.Hs.eg.db, keys = entrez, columns = "ALIAS", keytype = "ENTREZID")$ALIAS)
      match <- intersect(aliases, sc_universe)
      if (length(match) > 0) return(list(matched_gene = match[1], type = "OrgDb_Alias"))
    }
  }, error = function(e) NULL)
  
  # 5. BioMart (Database remoto esteso - NUOVO)
  # Usiamo la lista pre-fetched 'biomart_lookup'
  if (input_gene %in% names(biomart_lookup)) {
    bm_candidates <- biomart_lookup[[input_gene]]
    match <- intersect(bm_candidates, sc_universe)
    if (length(match) > 0) return(list(matched_gene = match[1], type = "BioMart_Synonym"))
  }
  
  return(list(matched_gene = NA, type = "Not_Found"))
}

process_gene_list <- function(gene_list, sc_data, universe) {
  results_list <- list()
  
  for (g in gene_list) {
    match_res <- find_gene_match(g, universe)
    
    if (!is.na(match_res$matched_gene)) {
      # FOUND: Filtra usando dplyr::filter esplicito
      rows <- sc_data %>% dplyr::filter(gene == match_res$matched_gene)
      
      rows <- rows %>% mutate(
        node = g,
        gene_rds = gene, 
        match_type = match_res$type
      )
      results_list[[g]] <- rows
      
    } else {
      # MISSING
      missing_row <- tibble(
        node = g,
        gene_rds = NA,
        match_type = "Not_Found",
        cluster = NA,
        p_val = NA,
        avg_log2fc = NA
      )
      results_list[[g]] <- missing_row
    }
  }
  
  dplyr::bind_rows(results_list)
}

# 8. ESECUZIONE ----------------------------------------------------------------
cat("\nAnalisi Visit 1...\n")
res_visit1 <- process_gene_list(visit1_genes, sc_de_data, sc_genes_universe)

cat("Analisi Visit 2...\n")
res_visit2 <- process_gene_list(visit2_genes, sc_de_data, sc_genes_universe)

# 9. OUTPUT EXCEL --------------------------------------------------------------
# Colonne desiderate
desired_cols <- c("node", "gene_rds", "match_type", "cluster", "p_val", "avg_log2fc", "p_val_adj")

safe_select <- function(df) {
  cols_present <- intersect(desired_cols, colnames(df))
  df %>% dplyr::select(all_of(cols_present), everything())
}

# Split Found/Missing
v1_found <- res_visit1 %>% dplyr::filter(match_type != "Not_Found") %>% safe_select()
v1_missing <- res_visit1 %>% dplyr::filter(match_type == "Not_Found") %>% dplyr::select(node, match_type)

v2_found <- res_visit2 %>% dplyr::filter(match_type != "Not_Found") %>% safe_select()
v2_missing <- res_visit2 %>% dplyr::filter(match_type == "Not_Found") %>% dplyr::select(node, match_type)

# Salvataggio
cat("\nGenerazione Excel...\n")

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "Bottom")

addWorksheet(wb, "VISIT1_Found")
writeData(wb, "VISIT1_Found", v1_found, headerStyle = hs)

addWorksheet(wb, "VISIT1_Missing")
writeData(wb, "VISIT1_Missing", v1_missing, headerStyle = hs)

addWorksheet(wb, "VISIT2_Found")
writeData(wb, "VISIT2_Found", v2_found, headerStyle = hs)

addWorksheet(wb, "VISIT2_Missing")
writeData(wb, "VISIT2_Missing", v2_missing, headerStyle = hs)

saveWorkbook(wb, output_file_path, overwrite = TRUE)

cat("----------------------------------------------------------\n")
cat("SUCCESS! File salvato in:\n", output_file_path, "\n")
cat("----------------------------------------------------------\n")
