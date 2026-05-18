# ==============================================================================
# SCRIPT ANALISI LONGITUDINALE WGCNA - TOP 10% DRIVERS (FACET PER ORIGINE)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(readxl)

options(stringsAsFactors = FALSE)

# 1. Configurazione Percorsi (INVARIANTI)
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

# 2. Caricamento Metadati (INVARIANTI)
# ==============================================================================
if (!file.exists(metadata_file)) stop("File metadati mancante")
metadata <- read_excel(metadata_file)
meta_clean <- metadata %>%
  select(`Proteomics ID`, condition) %>%
  rename(SampleID_Clean = `Proteomics ID`, Status = condition) %>%
  mutate(SampleID_Clean = trimws(SampleID_Clean))

# ------------------------------------------------------------------------------
# FUNZIONE: Estrazione Top 10% Geni (Highest |GS|) - (INVARIANTI)
# ------------------------------------------------------------------------------
get_top_GS_genes <- function(result_path, visit_name, gs_percentile = 0.9) {
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
  gs_col_name <- paste0("GS.", trait_col_wgcna)
  
  # 2. Seleziona i geni di quel modulo e calcola |GS|
  module_genes <- geneInfo %>% 
    filter(moduleColor == top_module) %>%
    mutate(
      Abs_GS = abs(get(gs_col_name))
    )
  
  # 3. Calcola il 90° percentile di Abs_GS (quindi il Top 10%)
  gs_threshold <- quantile(module_genes$Abs_GS, probs = gs_percentile, na.rm = TRUE)
  
  # 4. Seleziona i geni che superano la soglia
  top_genes <- module_genes %>%
    filter(Abs_GS >= gs_threshold) %>%
    select(1, moduleColor, Abs_GS) %>% # Colonna 1 è solitamente GeneSymbol o ID
    rename(Gene = 1, Module = moduleColor) %>%
    mutate(
      Driver_Origin = visit_name,
      GS_Threshold = gs_threshold
    )
  
  return(top_genes)
}

# Funzione Estrazione Dati (INVARIANTI)
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

# 3. Esecuzione: Selezione Top Drivers (INVARIANTI)
# ==============================================================================

# A. Trova i Top 10% Geni per ogni visita
top_gs_defs <- bind_rows(lapply(visits, function(v) get_top_GS_genes(visit_paths[[v]], v)))

# B. Stampa i risultati (per controllo)
cat("\n==========================================================\n")
cat("    TOP 10% PROTEINS SELECTED (Highest |GS|)\n")
cat("==========================================================\n")
print(top_gs_defs %>% select(Driver_Origin, Module, Gene, Abs_GS, GS_Threshold) %>% arrange(Driver_Origin, desc(Abs_GS)))
cat("\nNumero totale di proteine selezionate:", nrow(top_gs_defs), "\n")
cat("==========================================================\n\n")

# C. Preparazione pool di geni unici per estrazione dati
unique_genes_pool <- unique(top_gs_defs$Gene)

# 4. Estrazione Dati e Calcolo Score
# ==============================================================================

# A. Estrazione Raw Data
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

# C. Mapping
data_mapped <- data_normalized %>%
  inner_join(top_gs_defs %>% select(Gene, Driver_Origin, Module), by = "Gene")  # Aggiungiamo anche il modulo

# D. Aggregazione per Plot (Media di Gruppo per ogni singola proteina)
# L'aggregazione avviene solo per i gruppi SampleID_Clean (campioni),
# mantenendo Gene, Visit, Status e Driver_Origin distinti.
plot_data_aggregated <- data_mapped %>%
  group_by(Gene, Module, Driver_Origin, Visit, Status) %>%
  summarise(
    Mean = mean(Z_Score, na.rm = TRUE),
    SE = sd(Z_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Ordinamento fattori per coerenza visiva
plot_data_aggregated$Visit <- factor(plot_data_aggregated$Visit, levels = visits)
plot_data_aggregated$Driver_Origin <- factor(plot_data_aggregated$Driver_Origin, levels = visits)
plot_data_aggregated$Status <- factor(plot_data_aggregated$Status, levels = c("healed", "not_healed"))


# 5. Plotting (FACET PER DRIVER ORIGIN)
# ==============================================================================
# 5. Plotting (FACET PER DRIVER ORIGIN)
# ==============================================================================

# Generiamo un nome univoco per ogni proteina-modulo per distinguerli nel plot (Già Fatto)
# plot_data_aggregated <- plot_data_aggregated %>%
# mutate(Gene_Module = paste0(Gene, " (", Module, ")"))

p_facet <- ggplot(plot_data_aggregated, aes(x = Visit, y = Mean)) +
  
  # A. Facet: Separa i grafici in colonne in base alla visita da cui i driver sono stati selezionati
  facet_wrap(~ Driver_Origin, ncol = 3) +
  
  # B. Linee: Una linea per ogni combinazione Gene_Module e Status.
  # Grouping: interaction(Status, Gene_Module) assicura che ci sia una linea per (healed, CCL18) e una per (not_healed, CCL18)
  geom_line(aes(color = Status, linetype = Status, group = interaction(Status, Gene_Module)), linewidth = 0.8) +
  
  # C. Punti: Usiamo lo stesso raggruppamento per i punti
  geom_point(aes(color = Status, group = interaction(Status, Gene_Module)), size = 2, alpha = 0.7) +
  
  # D. Barre di errore (SE)
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE, color = Status, group = interaction(Status, Gene_Module)), width = 0.1, alpha = 0.3) +
  
  # E. Colori (Status)
  scale_color_manual(values = c("healed" = "#E41A1C", "not_healed" = "#377EB8")) +
  
  # F. Tipo di linea (Status)
  scale_linetype_manual(values = c("healed" = "solid", "not_healed" = "dashed")) +
  
  # G. Rendi le singole proteine leggermente trasparenti per evitare confusione
  # e aggiungi una leggera etichetta per identificare le linee (Opzionale, ma utile)
  geom_line(aes(group = Gene_Module), color = "gray50", linewidth = 0.1, alpha = 0.3) + # Linea di fondo grigia per ogni proteina
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Temporal analysis - Top Driver 10% |GS|",
    y = "Z-Score",
    x = "Visit"
  )

print(p_facet)

# Salvataggio Plot
dir.create(file.path(metadata_dir, "temporal_analysis"), showWarnings = FALSE)
ggsave(file.path(metadata_dir, "temporal_analysis", "temporal_analysis_top_10_percent_GS_Facet_Origin.pdf"), plot=p_facet, width=15, height=8)