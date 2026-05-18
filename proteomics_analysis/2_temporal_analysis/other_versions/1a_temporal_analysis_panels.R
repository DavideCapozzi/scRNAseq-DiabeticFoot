# ==============================================================================
# SCRIPT OTTIMALE: Z-SCORE (per Pannelli Diversi) + FACETING (per Chiarezza)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)

options(stringsAsFactors = FALSE)

# 1. Configurazione Percorsi
# ==============================================================================
base_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"
metadata_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined"
metadata_file <- file.path(metadata_dir, "condition_sample_match.xlsx")

visits <- c("V1", "V2", "V3")
visit_paths <- list(
  "V1" = file.path(base_dir, "DiabeticFoot_V1", "Results_V1"),
  "V2" = file.path(base_dir, "DiabeticFoot_V2", "Results_V2"),
  "V3" = file.path(base_dir, "DiabeticFoot_V3", "Results_V3")
)

trait_col_wgcna <- "healed" 

# 2. Caricamento Metadati
# ==============================================================================
if (!file.exists(metadata_file)) stop(paste("File metadati non trovato:", metadata_file))
metadata <- read_excel(metadata_file)

if (!("Proteomics ID" %in% names(metadata)) | !("condition" %in% names(metadata))) {
  stop("Colonne 'Proteomics ID' e 'condition' mancanti nell'Excel.")
}

meta_clean <- metadata %>%
  select(`Proteomics ID`, condition) %>%
  rename(SampleID_Clean = `Proteomics ID`, Status = condition) %>%
  mutate(SampleID_Clean = trimws(SampleID_Clean))

# 3. Funzioni Helper
# ==============================================================================

get_top_hub_genes <- function(result_path, visit_name) {
  path1 <- file.path(result_path, "txtFile", "geneInfo.txt")
  path2 <- file.path(result_path, "geneInfo.txt")
  
  if (file.exists(path1)) gene_info_file <- path1
  else if (file.exists(path2)) gene_info_file <- path2
  else stop(paste("geneInfo.txt non trovato in:", result_path))
  
  geneInfo <- read.table(gene_info_file, header = TRUE, sep = "\t")
  
  module_stats <- geneInfo %>%
    filter(moduleColor != "grey") %>%
    group_by(moduleColor) %>%
    summarise(mean_GS = mean(abs(get(paste0("GS.", trait_col_wgcna))), na.rm = TRUE)) %>%
    arrange(desc(mean_GS))
  
  if(nrow(module_stats) == 0) return(NULL)
  top_module <- module_stats$moduleColor[1]
  
  mm_col <- paste0("MM.", top_module)
  gs_col <- paste0("GS.", trait_col_wgcna)
  
  top_genes <- geneInfo %>%
    filter(moduleColor == top_module) %>%
    mutate(HubScore = abs(get(gs_col)) + abs(get(mm_col))) %>%
    arrange(desc(HubScore)) %>%
    head(2) %>%
    pull(1) 
  
  return(data.frame(Gene = top_genes, Driver_Origin = visit_name, Module = top_module))
}

get_expression_data <- function(result_path, visit_name, target_genes, metadata_df) {
  rdata_1 <- file.path(result_path, "RData", "dataInput.RData")
  rdata_2 <- file.path(result_path, "dataInput.RData")
  
  if (file.exists(rdata_1)) rdata_file <- rdata_1
  else if (file.exists(rdata_2)) rdata_file <- rdata_2
  else stop(paste("dataInput.RData non trovato per", visit_name))
  
  env <- new.env()
  load(rdata_file, envir = env)
  datExpr <- env$datExpr
  
  common_genes <- intersect(colnames(datExpr), target_genes)
  if (length(common_genes) == 0) return(NULL)
  
  expr_subset <- datExpr[, common_genes, drop = FALSE]
  
  df_long <- as.data.frame(expr_subset) %>%
    rownames_to_column("SampleID_Clean") %>%
    mutate(SampleID_Clean = trimws(SampleID_Clean)) %>% 
    pivot_longer(cols = all_of(common_genes), names_to = "Gene", values_to = "Expression") %>%
    mutate(Visit = visit_name)
  
  df_merged <- df_long %>%
    left_join(metadata_df, by = "SampleID_Clean") %>%
    filter(!is.na(Status))
  
  return(df_merged)
}

# 4. Esecuzione
# ==============================================================================

# A. Identifica geni
all_target_genes_info <- data.frame()
for (v in visits) {
  tryCatch({
    top_v <- get_top_hub_genes(visit_paths[[v]], v)
    if(!is.null(top_v)) all_target_genes_info <- rbind(all_target_genes_info, top_v)
  }, error = function(e) message(paste("Skip ricerca", v, ":", e$message)))
}

if(nrow(all_target_genes_info) == 0) stop("Nessun gene trovato.")
unique_genes <- unique(all_target_genes_info$Gene)

# B. Estrai dati (Grezzi per ora)
longitudinal_data <- data.frame()
for (v in visits) {
  tryCatch({
    curr_data <- get_expression_data(visit_paths[[v]], v, unique_genes, meta_clean)
    if(!is.null(curr_data)) longitudinal_data <- rbind(longitudinal_data, curr_data)
  }, error = function(e) message(paste("Skip dati", v, ":", e$message)))
}

# C. Calcolo Z-Score (IMPORTANTE per pannelli diversi)
# Normalizziamo ogni gene indipendentemente. In questo modo portiamo
# proteine con range 0-10 e proteine con range 1000-5000 sulla stessa scala relativa.
plot_data <- longitudinal_data %>%
  group_by(Gene) %>%
  mutate(Z_Score = scale(Expression)) %>% 
  ungroup()

# D. Aggregazione Media
gene_origins <- all_target_genes_info %>%
  group_by(Gene) %>%
  slice(1) %>%
  mutate(Label = paste0(Gene, "\n(Driver: ", Driver_Origin, ")")) %>%
  select(Gene, Label)

plot_data_summary <- plot_data %>%
  left_join(gene_origins, by = "Gene") %>%
  group_by(Gene, Visit, Status, Label) %>%
  summarise(
    Mean_Z = mean(Z_Score, na.rm = TRUE),
    SE_Z = sd(Z_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

plot_data_summary$Visit <- factor(plot_data_summary$Visit, levels = c("V1", "V2", "V3"))

# 5. Plotting (Z-Score + Faceting)
# ==============================================================================

p <- ggplot(plot_data_summary, aes(x = Visit, y = Mean_Z, group = Status, color = Status)) +
  
  # Linee e Punti
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  # Barre di errore
  geom_errorbar(aes(ymin = Mean_Z - SE_Z, ymax = Mean_Z + SE_Z), 
                width = 0.1, alpha = 0.7) +
  
  # PANNELLI: Usiamo facet_wrap. 
  # Qui 'scales="fixed"' è appropriato perché sono Z-score (tutti sulla stessa scala),
  # rendendo i pannelli perfettamente confrontabili tra loro.
  facet_wrap(~Label, ncol = 3) +
  
  # Colori "Traffic Light" (Verde=Guarito, Rosso=Non Guarito)
  scale_color_manual(values = c("healed" = "#009E73", "not_healed" = "#D55E00"), 
                     labels = c("Healed", "Not Healed")) +
  
  theme_bw() +
  
  # Linea dello zero per riferimento (media della popolazione)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  labs(
    title = "Profili Temporali Normalizzati (Z-Score)",
    y = "Z-Score Expression (Mean +/- SE)",
    x = "Time Point",
    color = "Clinical Status"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.background = element_rect(fill = "#e6e6e6"), 
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(p)

# # Salvataggio
# output_file <- file.path(metadata_dir, "Trajectory_ZScore_Final.pdf")
# ggsave(output_file, plot = p, width = 12, height = 10)
# print(paste("Grafico salvato in:", output_file))