rm(list = ls())
library(enrichR)
library(ggplot2)
library(xlsx)

# ---------------------------------------------------------------
# PARAMETRI GLOBALI
# ---------------------------------------------------------------
pval_thr    <- 0.1
# Definisco i DB: includo sia DisGeNET che GO per evitare errori di estrazione
dbs         <- c("DisGeNET", "GO_Biological_Process_2023")
all_modules <- c("turquoise", "blue", "brown", "green", "yellow", "red", "black")
visits      <- c("V1", "V2", "V3")

# Cartella di output
res_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/enrichment_res/DisGeNET"

if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
} else { print("the dir already exists")
 }

# Connessione a Enrichr
setEnrichrSite("Enrichr")

# ---------------------------------------------------------------
# LOOP ESTERNO: visite
# ---------------------------------------------------------------
for (visit in visits) {
  
  cat("\n\n##############################################\n")
  cat("VISIT:", visit, "\n")
  cat("##############################################\n")
  
  # Path input aggiornato
  input_path <- paste0(
    "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/",
    "sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/",
    "CeciliaMorgantini/dataset/DiabeticFoot_", visit,
    "/Results_", visit, "/txtFile"
  )
  
  visit_res_path <- file.path(res_path, visit)
  if (!dir.exists(visit_res_path)) dir.create(visit_res_path, recursive = TRUE)
  
  gene_file <- file.path(input_path, "geneInfo.txt")
  if (!file.exists(gene_file)) {
    cat("  -> FILE NOT FOUND:", gene_file, "- Skipping visit", visit, "\n")
    next
  }
  
  df <- read.delim(gene_file, header = TRUE, sep = "\t")
  
  # ---------------------------------------------------------------
  # LOOP INTERNO: moduli
  # ---------------------------------------------------------------
  for (mod in all_modules) {
    
    cat("\n  Module:", mod, "| Visit:", visit, "\n")
    
    # Selezione geni (usando la colonna corretta del tuo geneInfo.txt)
    # Assumo che la prima colonna siano i simboli e 'moduleColor' sia la colonna dei colori
    genes_to_test <- df$geneSymbol[df$moduleColor == mod] # Verifica se la colonna si chiama geneSymbol o gene
    if(is.null(genes_to_test)) genes_to_test <- df[df$moduleColor == mod, 1] 
    
    genes_to_test <- na.omit(genes_to_test)
    
    if (length(genes_to_test) < 3) { # Richiesti almeno 3 geni per un senso biologico
      cat("  -> Too few genes - Skipping.\n")
      next
    }
    
    # Enrichment
    enriched <- tryCatch(
      enrichr(genes_to_test, dbs),
      error = function(e) return(NULL)
    )
    
    if (is.null(enriched)) {
      cat("  -> Enrichr connection error - Skipping.\n")
      next
    }
    
    # --- PROCESSO DISGENET ---
    res_dis <- enriched[["DisGeNET"]]
    res_go <- res_dis
    
    # --- PROCESSO GO (Per i Plot) ---
    #res_go <- enriched[["GO_Biological_Process_2023"]]
    
    # Funzione di filtraggio e validazione robusta
    process_res <- function(res_table) {
      if (is.null(res_table) || !is.data.frame(res_table) || nrow(res_table) == 0) {
        return(NULL)
      }
      filtered <- res_table[res_table$Adjusted.P.value < pval_thr, ]
      if (nrow(filtered) == 0) return(NULL)
      return(filtered)
    }
    
    res_filtered_go <- process_res(res_go)
    
    # Controllo di sicurezza per evitare l'errore "argument is of length zero"
    # if (is.null(res_filtered_go)) {
    #   cat("  -> No significant terms found for this module.\n")
    #   next
    # }

    #cat("  -> Significant terms found:", nrow(res_filtered_go), "\n")

    # 1. Salvataggio Excel (DisGeNET se presente, altrimenti GO)
    res_to_save <- if(!is.null(process_res(res_dis))) process_res(res_dis) else res_filtered_go
    xlsx_file <- file.path(visit_res_path, paste0(visit, "_", mod, "_DisGeNET_enrichment.xlsx"))
    write.xlsx(res_to_save, xlsx_file, row.names = FALSE)
    
    # 2. Preparazione Top 10 per Plot
    res_top10 <- head(res_filtered_go[order(res_filtered_go$Adjusted.P.value), ], 10)
    res_top10$Count <- as.numeric(sub("/.*", "", res_top10$Overlap))
    
    # 3. Plot
    p <- ggplot(res_top10, aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
      geom_point(aes(size = Count, color = Adjusted.P.value)) +
      scale_color_gradient(low = "red", high = "blue") +
      scale_size_continuous(range = c(5, 12)) +
      theme_classic() +
      labs(
        title = paste0("Top 10 DisGeNET terms (", tools::toTitleCase(mod), " - ", visit, ")"),
        subtitle = paste0("Adj. p-value < ", pval_thr),
        x = "Combined Score", y = ""
      )
    
    # 4. Salvataggio PDF
    pdf_file <- file.path(visit_res_path, paste0(visit, "_", mod, "_plot.pdf"))
    ggsave(filename = pdf_file, plot = p, width = 10, height = 7)
  }
}
cat("\nProcesso completato.\n")
