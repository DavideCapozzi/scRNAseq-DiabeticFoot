# ==============================================================================
# SCRIPT ANALISI LONGITUDINALE WGCNA - TOP 3 DRIVERS (GS + MM)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(readxl)

options(stringsAsFactors = FALSE)

# 1. Configurazione Percorsi
# ==============================================================================
metadata_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined"
dataset_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"
metadata_file <- file.path(metadata_dir, "condition_sample_match.xlsx")

visits <- c("V1", "V2", "V3")
visit_paths <- list(
  "V1" = file.path(dataset_dir, "DiabeticFoot_V1", "Results_V1"),
  "V2" = file.path(dataset_dir, "DiabeticFoot_V2", "Results_V2"),
  "V3" = file.path(dataset_dir, "DiabeticFoot_V3", "Results_V3")
)
trait_col_wgcna <- "healed" 

# 2. Caricamento Metadati
# ==============================================================================
if (!file.exists(metadata_file)) stop("File metadati mancante")
metadata <- read_excel(metadata_file)
meta_clean <- metadata %>%
  select(`Proteomics ID`, condition) %>%
  rename(SampleID_Clean = `Proteomics ID`, Status = condition) %>%
  mutate(SampleID_Clean = trimws(SampleID_Clean))

# ------------------------------------------------------------------------------
# NUOVA FUNZIONE: Estrazione Top 3 Geni (Highest GS + MM)
# ------------------------------------------------------------------------------
get_top_hub_genes <- function(result_path, visit_name, top_n = 3) {
  path1 <- file.path(result_path, "txtFile", "geneInfo.txt")
  path2 <- file.path(result_path, "geneInfo.txt")
  if (file.exists(path1)) f <- path1 else if (file.exists(path2)) f <- path2 else return(NULL)
  
  geneInfo <- read.table(f, header = TRUE, sep = "\t")
  
  # 1. Identifica il modulo più correlato al trait (Top Module)
  module_stats <- geneInfo %>%
    filter(moduleColor != "grey") %>%
    group_by(moduleColor) %>%
    summarise(mean_GS = mean(abs(get(paste0("GS.", trait_col_wgcna))), na.rm = TRUE)) %>%
    arrange(desc(mean_GS))
  
  if(nrow(module_stats) == 0) return(NULL)
  top_module <- module_stats$moduleColor[1]
  
  # 2. Seleziona i geni di quel modulo
  module_genes <- geneInfo %>% 
    filter(moduleColor == top_module)
  
  # 3. Calcola Hub Score = |GS| + |MM|
  # Serve costruire dinamicamente il nome della colonna MM (es. MM.blue)
  mm_col_name <- paste0("MM.", top_module)
  gs_col_name <- paste0("GS.", trait_col_wgcna)
  
  top_genes <- module_genes %>%
    mutate(
      Abs_GS = abs(get(gs_col_name)),
      Abs_MM = abs(get(mm_col_name)),
      Hub_Score = Abs_GS + Abs_MM
    ) %>%
    arrange(desc(Hub_Score)) %>% # Ordina dal più alto al più basso
    slice_head(n = top_n) %>%    # Prendi i primi 3
    select(1, moduleColor, Abs_GS, Abs_MM, Hub_Score) %>% # Colonna 1 è solitamente GeneSymbol o ID
    rename(Gene = 1, Module = moduleColor) %>%
    mutate(Driver_Origin = visit_name)
  
  return(top_genes)
}

# Funzione Estrazione Dati (Invariata nella logica, usata dopo)
get_expression_data <- function(result_path, visit_name, target_genes, metadata_df) {
  rdata_1 <- file.path(result_path, "RData", "dataInput.RData")
  rdata_2 <- file.path(result_path, "dataInput.RData")
  if (file.exists(rdata_1)) f <- rdata_1 else if (file.exists(rdata_2)) f <- rdata_2 else return(NULL)
  
  env <- new.env(); load(f, envir = env); datExpr <- env$datExpr
  common <- intersect(colnames(datExpr), target_genes)
  if (length(common) == 0) return(NULL)
  
  df <- as.data.frame(datExpr[, common, drop=F]) %>%
    rownames_to_column("SampleID_Clean") %>%
    mutate(SampleID_Clean = trimws(SampleID_Clean)) %>%
    pivot_longer(-SampleID_Clean, names_to="Gene", values_to="Expression") %>%
    mutate(Visit = visit_name) %>%
    left_join(metadata_df, by="SampleID_Clean") %>%
    filter(!is.na(Status))
  return(df)
}

# 3. Esecuzione: Selezione Top Drivers
# ==============================================================================

# A. Trova i Top 3 Geni per ogni visita
top3_defs <- bind_rows(lapply(visits, function(v) get_top_hub_genes(visit_paths[[v]], v)))

# B. Stampa i risultati (COME RICHIESTO)
cat("\n==========================================================\n")
cat("   TOP 3 PROTEINS SELECTED (Highest GS + MM)\n")
cat("==========================================================\n")
print(top3_defs %>% select(Driver_Origin, Module, Gene, Abs_GS, Abs_MM, Hub_Score))
cat("\n==========================================================\n\n")

# C. Preparazione pool di geni unici per estrazione dati
# (Nota: se un gene è top in V1 e V2, qui appare una volta sola per l'estrazione raw data)
unique_genes_pool <- unique(top3_defs$Gene)

# 4. Estrazione Dati e Calcolo Score
# ==============================================================================

# A. Estrazione Raw Data (Dati di espressione per tutte le visite)
raw_data <- data.frame()
for (v in visits) {
  tryCatch({
    curr_data <- get_expression_data(visit_paths[[v]], v, unique_genes_pool, meta_clean)
    if(!is.null(curr_data)) raw_data <- rbind(raw_data, curr_data)
  }, error = function(e) message(paste("Err data", v, ":", e$message)))
}

# B. Normalizzazione (Z-Score)
data_normalized <- raw_data %>%
  group_by(Gene) %>%
  mutate(Z_Score = as.numeric(scale(Expression))) %>%
  ungroup()

# C. Mapping e Aggregazione (Punto cruciale)
# Qui facciamo il JOIN con la tabella 'top3_defs'.
# Se il gene "X" è presente sia in "V1" che in "V2" in top3_defs,
# il join duplicherà le righe di espressione: una assegnata a Driver_Origin V1, una a V2.
data_mapped <- data_normalized %>%
  inner_join(top3_defs %>% select(Gene, Driver_Origin), by = "Gene") 

# D. Aggregazione per Paziente
patient_hub_scores <- data_mapped %>%
  group_by(SampleID_Clean, Visit, Status, Driver_Origin) %>%
  summarise(
    Patient_Module_Score = mean(Z_Score, na.rm = TRUE),
    .groups = 'drop'
  )

# E. Aggregazione Finale per il Plot (Media di Gruppo)
plot_data <- patient_hub_scores %>%
  group_by(Visit, Status, Driver_Origin) %>%
  summarise(
    Mean = mean(Patient_Module_Score, na.rm = TRUE),
    SE = sd(Patient_Module_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Ordinamento fattori per coerenza visiva
plot_data$Visit <- factor(plot_data$Visit, levels = c("V1", "V2", "V3"))
plot_data$Driver_Origin <- factor(plot_data$Driver_Origin, levels = c("V1", "V2", "V3"))
plot_data$Status <- factor(plot_data$Status, levels = c("healed", "not_healed"))

# 5. Plotting
# ==============================================================================

p <- ggplot(plot_data, aes(x = Visit, y = Mean, group = interaction(Driver_Origin, Status))) +
  
  # A. Barre di errore
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE, color = Driver_Origin), 
                width = 0.1, alpha = 0.3) +
  
  # B. Linee
  geom_line(aes(color = Driver_Origin, alpha = Status, linewidth = Status)) +
  
  # C. Punti
  geom_point(aes(color = Driver_Origin, alpha = Status), size = 3) +
  
  # D. Colori (Origin)
  scale_color_manual(values = c("V1" = "red", "V2" = "blue", "V3" = "green"),
                     name = "Visit") +
  
  # E. Trasparenza (Status)
  scale_alpha_manual(values = c("healed" = 1, "not_healed" = 0.4),
                     name = "Condition") +
  
  # F. Spessore (Status)
  scale_linewidth_manual(values = c("healed" = 1.5, "not_healed" = 0.8),
                         name = "Condition") +
  
  theme_bw() +
  labs(
    title = "Temporal Analysis - top 3 driver proteins",
    subtitle = "Drivers selected per visit based on GS and MM",
    y = "Average Z-Score",
    x = "Visit"
  )

print(p)

# Salvataggio Plot
ggsave(file.path(metadata_dir, "temporal_analysis", "temporal_analysis_top_3_GS_MM.pdf"), plot=p, width=10, height=8)

# # Salvataggio Lista Geni
# write.csv(top3_defs, file.path(metadata_dir, "Top3_Hub_Genes_List.csv"), row.names = FALSE)