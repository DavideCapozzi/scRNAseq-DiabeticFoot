# ==============================================================================
# SCRIPT ANALISI LONGITUDINALE WGCNA (Colori Specifici per Gene)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(scales) # Per la gestione automatica dei colori

options(stringsAsFactors = FALSE)

# 1. Configurazione Percorsi
# ==============================================================================

# Directory contenente il file di metadati (analysis_joined)
metadata_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined"

# Directory radice dei dataset WGCNA
dataset_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"

# Percorso completo file metadati
metadata_file <- file.path(metadata_dir, "condition_sample_match.xlsx")

# Definizione visite
visits <- c("V1", "V2", "V3")
visit_paths <- list(
  "V1" = file.path(dataset_dir, "DiabeticFoot_V1", "Results_V1"),
  "V2" = file.path(dataset_dir, "DiabeticFoot_V2", "Results_V2"),
  "V3" = file.path(dataset_dir, "DiabeticFoot_V3", "Results_V3")
)

trait_col_wgcna <- "healed" 

# 2. Caricamento Metadati
# ==============================================================================

if (!file.exists(metadata_file)) {
  stop(paste("File metadati non trovato:", metadata_file))
}

print("Caricamento metadati...")
metadata <- read_excel(metadata_file)

# Standardizzazione nomi colonne
if (!("Proteomics ID" %in% names(metadata)) | !("condition" %in% names(metadata))) {
  stop("Le colonne 'Proteomics ID' e 'condition' sono richieste nel file Excel.")
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
  
  # Trova modulo top
  module_stats <- geneInfo %>%
    filter(moduleColor != "grey") %>%
    group_by(moduleColor) %>%
    summarise(mean_GS = mean(abs(get(paste0("GS.", trait_col_wgcna))), na.rm = TRUE)) %>%
    arrange(desc(mean_GS))
  
  if(nrow(module_stats) == 0) return(NULL)
  top_module <- module_stats$moduleColor[1]
  
  # Trova top 3 geni
  mm_col <- paste0("MM.", top_module)
  gs_col <- paste0("GS.", trait_col_wgcna)
  
  top_genes <- geneInfo %>%
    filter(moduleColor == top_module) %>%
    mutate(HubScore = abs(get(gs_col)) + abs(get(mm_col))) %>%
    arrange(desc(HubScore)) %>%
    head(3) %>%
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
  
  # Join con metadati
  df_merged <- df_long %>%
    left_join(metadata_df, by = "SampleID_Clean") %>%
    filter(!is.na(Status))
  
  return(df_merged)
}

# 4. Esecuzione e Preparazione Dati
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

# B. Estrai dati
longitudinal_data <- data.frame()
for (v in visits) {
  tryCatch({
    curr_data <- get_expression_data(visit_paths[[v]], v, unique_genes, meta_clean)
    if(!is.null(curr_data)) longitudinal_data <- rbind(longitudinal_data, curr_data)
  }, error = function(e) message(paste("Skip dati", v, ":", e$message)))
}

# C. Calcolo Z-Score
plot_data <- longitudinal_data %>%
  group_by(Gene) %>%
  mutate(Z_Score = scale(Expression)) %>% 
  ungroup()

# D. Aggregazione Media
gene_origins <- all_target_genes_info %>%
  group_by(Gene) %>%
  slice(1) %>%
  select(Gene, Driver_Origin)

plot_data_summary <- plot_data %>%
  left_join(gene_origins, by = "Gene") %>%
  group_by(Gene, Visit, Status, Driver_Origin) %>%
  summarise(
    Mean_Z = mean(Z_Score, na.rm = TRUE),
    SE_Z = sd(Z_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

plot_data_summary$Visit <- factor(plot_data_summary$Visit, levels = c("V1", "V2", "V3"))

# Creiamo la variabile combinata per il mapping dei colori
plot_data_summary$Gene_Status <- paste(plot_data_summary$Gene, plot_data_summary$Status, sep = "_")

# 5. Generazione Palette Colori Dinamica (Coppie)
# ==============================================================================

# Otteniamo i geni unici ordinati per visita di origine
ordered_genes <- plot_data_summary %>%
  arrange(Driver_Origin, Gene) %>%
  pull(Gene) %>%
  unique()

n_genes <- length(ordered_genes)

# Genera n colori base ben distinti
base_hues <- hue_pal()(n_genes)

# Costruiamo il vettore dei colori manuali
custom_colors <- c()

for (i in 1:n_genes) {
  gene <- ordered_genes[i]
  base_col <- base_hues[i]
  
  # Definiamo i nomi delle chiavi
  key_healed <- paste(gene, "healed", sep="_")
  key_not_healed <- paste(gene, "not_healed", sep="_")
  
  # Assegnazione colori:
  # Healed = Colore Base Pieno
  # Not Healed = Colore Base Schiarito (più pastello)
  
  # Funzione semplice per schiarire il colore in esadecimale
  light_col <- grDevices::adjustcolor(base_col, alpha.f = 0.4) 
  
  custom_colors[key_healed] <- base_col
  custom_colors[key_not_healed] <- light_col
}

# 6. Plotting
# ==============================================================================

p <- ggplot(plot_data_summary, aes(x = Visit, y = Mean_Z, group = Gene_Status)) +
  
  # Linee e Punti
  geom_line(aes(color = Gene_Status, linetype = Driver_Origin), linewidth = 1) +
  geom_point(aes(color = Gene_Status, shape = Driver_Origin), size = 3) +
  
  # Barre di errore
  geom_errorbar(aes(ymin = Mean_Z - SE_Z, ymax = Mean_Z + SE_Z, color = Gene_Status), 
                width = 0.1, alpha = 0.6) +
  
  # Applicazione Palette Dinamica
  scale_color_manual(values = custom_colors) +
  
  # Stili linee
  scale_linetype_manual(values = c("V1" = "solid", "V2" = "longdash", "V3" = "dotted"),
                        name = "Driver Origin") +
  scale_shape_manual(values = c("V1" = 16, "V2" = 17, "V3" = 15),
                     name = "Driver Origin") +
  
  theme_bw() +
  
  labs(
    title = "Traiettoria Longitudinale Geni Driver (Per Proteina)",
    subtitle = "Ogni colore rappresenta una proteina (Pieno=Healed, Chiaro=Not Healed)",
    y = "Z-score Expression",
    x = "Time Point",
    color = "Protein & Status"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right",
    # Riduciamo la dimensione del testo della legenda perché ci saranno molte voci
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.4, "cm")
  ) +
  # Organizziamo la legenda colori su due colonne per risparmiare spazio verticale
  guides(color = guide_legend(ncol = 2))

print(p)

# # Salvataggio
# output_file <- file.path(metadata_dir, "Trajectory_PairedColors.pdf")
# ggsave(output_file, plot = p, width = 14, height = 9) # Leggermente più largo per la legenda
# print(paste("Grafico salvato in:", output_file))