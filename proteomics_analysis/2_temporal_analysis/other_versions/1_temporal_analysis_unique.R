# ==============================================================================
# SCRIPT ANALISI LONGITUDINALE WGCNA - 6 LINEE (Colori per Origine, Stile per Status)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(readxl)

options(stringsAsFactors = FALSE)

# 1. Configurazione Percorsi (Invariata)
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

# 2. Caricamento Metadati & Funzioni (Invariata)
# ==============================================================================
# (Eseguire questa parte come nello script precedente per caricare metadata e definire funzioni)

if (!file.exists(metadata_file)) stop("File metadati mancante")
metadata <- read_excel(metadata_file)
meta_clean <- metadata %>%
  select(`Proteomics ID`, condition) %>%
  rename(SampleID_Clean = `Proteomics ID`, Status = condition) %>%
  mutate(SampleID_Clean = trimws(SampleID_Clean))

# Funzione Estrazione Geni
get_module_genes <- function(result_path, visit_name) {
  path1 <- file.path(result_path, "txtFile", "geneInfo.txt")
  path2 <- file.path(result_path, "geneInfo.txt")
  if (file.exists(path1)) f <- path1 else if (file.exists(path2)) f <- path2 else return(NULL)
  
  geneInfo <- read.table(f, header = TRUE, sep = "\t")
  module_stats <- geneInfo %>%
    filter(moduleColor != "grey") %>%
    group_by(moduleColor) %>%
    summarise(mean_GS = mean(abs(get(paste0("GS.", trait_col_wgcna))), na.rm = TRUE)) %>%
    arrange(desc(mean_GS))
  
  if(nrow(module_stats) == 0) return(NULL)
  top_module <- module_stats$moduleColor[1]
  
  # Estrae la colonna 1 (GeneName)
  genes <- geneInfo %>% filter(moduleColor == top_module) %>% pull(1)
  return(data.frame(Gene = genes, Driver_Origin = visit_name, Module = top_module))
}

# Funzione Estrazione Dati
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

# 3. Preparazione Dati
# ==============================================================================
# A. Definizioni Moduli
all_module_defs <- bind_rows(lapply(visits, function(v) get_module_genes(visit_paths[[v]], v)))
unique_genes_pool <- unique(all_module_defs$Gene)

# B. Dati Raw
raw_data <- bind_rows(lapply(visits, function(v) get_expression_data(visit_paths[[v]], v, unique_genes_pool, meta_clean)))

# C. Aggregazione
plot_data <- raw_data %>%
  group_by(Gene) %>%
  mutate(Z_Score = as.numeric(scale(Expression))) %>% # Z-score gene-wise
  ungroup() %>%
  inner_join(all_module_defs, by = "Gene") %>% # Duplica geni se presenti in più firme
  group_by(SampleID_Clean, Visit, Status, Driver_Origin) %>%
  summarise(Patient_Score = mean(Z_Score, na.rm=TRUE), .groups="drop") %>% # Media per paziente
  group_by(Visit, Status, Driver_Origin) %>%
  summarise(
    Mean = mean(Patient_Score, na.rm=TRUE),
    SE = sd(Patient_Score, na.rm=TRUE)/sqrt(n()),
    .groups="drop"
  )

# Ordinamento Fattori
plot_data$Visit <- factor(plot_data$Visit, levels = c("V1", "V2", "V3"))
plot_data$Driver_Origin <- factor(plot_data$Driver_Origin, levels = c("V1", "V2", "V3"))
plot_data$Status <- factor(plot_data$Status, levels = c("healed", "not_healed")) # Ordine legenda

# 4. Plotting (Colori per Origine, Linee per Status)
# ==============================================================================
# ==============================================================================
# 4. Elaborazione Dati con FILTRO ESCLUSIVITÀ
# ==============================================================================

# A. Identifica le liste di geni per ogni visita
all_module_defs <- data.frame()
for (v in visits) {
  tryCatch({
    mod_v <- get_module_genes(visit_paths[[v]], v)
    if(!is.null(mod_v)) all_module_defs <- rbind(all_module_defs, mod_v)
  }, error = function(e) message(paste("Err", v, ":", e$message)))
}

# --- NUOVO PASSAGGIO: FILTRO PER UNICITÀ ASSOLUTA ---
print("Conteggio geni pre-filtro:")
print(table(all_module_defs$Driver_Origin))

# Identifichiamo i geni che sono presenti in più di una visita
gene_counts <- all_module_defs %>%
  group_by(Gene) %>%
  summarise(Appearances = n())

# Teniamo SOLO i geni che appaiono esattamente 1 volta (Esclusivi)
exclusive_genes <- gene_counts %>%
  filter(Appearances == 1) %>%
  pull(Gene)

# Filtriamo il dataset delle definizioni
all_module_defs_exclusive <- all_module_defs %>%
  filter(Gene %in% exclusive_genes)


# ==============================================================================
# ANALISI SOVRAPPOSIZIONE GENI (VENN DIAGRAM)
# ==============================================================================

# 1. Preparazione della Lista per il Venn
# Creiamo una lista dove ogni elemento contiene i geni di una visita
venn_list <- list(
  V1 = all_module_defs %>% filter(Driver_Origin == "V1") %>% pull(Gene),
  V2 = all_module_defs %>% filter(Driver_Origin == "V2") %>% pull(Gene),
  V3 = all_module_defs %>% filter(Driver_Origin == "V3") %>% pull(Gene)
)

# 2. Creazione e Salvataggio del Diagramma di Venn
# Nota: Richiede il pacchetto ggvenn. Se non ce l'hai: install.packages("ggvenn")
if(require("ggvenn", quietly = TRUE)) {
  p_venn <- ggvenn(
    venn_list, 
    fill_color = c("green", "blue2", "turquoise"),
    stroke_size = 0.5, 
    show_percentage = F ,
    set_name_size = 4
  ) + ggtitle("Proteins overlap between visits")
  
  print(p_venn)
  # Opzionale: ggsave(file.path(metadata_dir, "Venn_Overlap.pdf"), p_venn)
} else {
  message("Pacchetto 'ggvenn' non trovato. Salto il plot del Venn.")
}

# 3. Calcolo Liste: Unici e Comuni
# Intersezione di TUTTI e 3 (Geni presenti in V1, V2 e V3)
common_all <- Reduce(intersect, venn_list)

# Geni Unici (Esclusivi) per ogni visita
unique_V1 <- setdiff(venn_list$V1, union(venn_list$V2, venn_list$V3))
unique_V2 <- setdiff(venn_list$V2, union(venn_list$V1, venn_list$V3))
unique_V3 <- setdiff(venn_list$V3, union(venn_list$V1, venn_list$V2))

# 4. Stampa dei Risultati in Console
cat("\n==============================================\n")
cat("      RISULTATI SOVRAPPOSIZIONE GENI\n")
cat("==============================================\n")

cat(paste0("\n--- GENI COMUNI A TUTTE LE VISITE (n=", length(common_all), ") ---\n"))
if(length(common_all) > 0) print(common_all) else cat("Nessun gene comune a tutte e 3 le visite.\n")

cat(paste0("\n--- GENI UNICI V1 (n=", length(unique_V1), ") ---\n"))
if(length(unique_V1) > 0) print(unique_V1) else cat("Nessun gene esclusivo per V1.\n")

cat(paste0("\n--- GENI UNICI V2 (n=", length(unique_V2), ") ---\n"))
if(length(unique_V2) > 0) print(unique_V2) else cat("Nessun gene esclusivo per V2.\n")

cat(paste0("\n--- GENI UNICI V3 (n=", length(unique_V3), ") ---\n"))
if(length(unique_V3) > 0) print(unique_V3) else cat("Nessun gene esclusivo per V3.\n")

cat("\n==============================================\n\n")

# ==============================================================================

print("--- FILTRO APPLICATO: RIMOSSI GENI CONDIVISI ---")
print("Conteggio geni post-filtro (Esclusivi per visita):")
print(table(all_module_defs_exclusive$Driver_Origin))

if(nrow(all_module_defs_exclusive) == 0) stop("Attenzione: Nessun gene è esclusivo! Tutti i geni sono condivisi tra le visite.")

# Aggiorniamo la pool di geni unici da estrarre
unique_genes_pool <- unique(all_module_defs_exclusive$Gene)

# B. Estrai i dati di espressione raw (solo per i geni esclusivi)
raw_data <- data.frame()
for (v in visits) {
  tryCatch({
    curr_data <- get_expression_data(visit_paths[[v]], v, unique_genes_pool, meta_clean)
    if(!is.null(curr_data)) raw_data <- rbind(raw_data, curr_data)
  }, error = function(e) message(paste("Err data", v, ":", e$message)))
}

# C. Calcolo Z-Score e Aggregazione
# --------------------------------------------------------------------------
# 1. Z-Score
data_normalized <- raw_data %>%
  group_by(Gene) %>%
  mutate(Z_Score = as.numeric(scale(Expression))) %>%
  ungroup()

# 2. JOIN con la lista ESCLUSIVA
data_mapped <- data_normalized %>%
  inner_join(all_module_defs_exclusive, by = "Gene") 

# 3. Aggregazione per Paziente
patient_module_scores <- data_mapped %>%
  group_by(SampleID_Clean, Visit, Status, Driver_Origin) %>%
  summarise(
    Patient_Module_Score = mean(Z_Score, na.rm = TRUE),
    .groups = 'drop'
  )

# 4. Aggregazione Finale per il Plot
plot_data <- patient_module_scores %>%
  group_by(Visit, Status, Driver_Origin) %>%
  summarise(
    Mean = mean(Patient_Module_Score, na.rm = TRUE),
    SE = sd(Patient_Module_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Ordinamento fattori
plot_data$Visit <- factor(plot_data$Visit, levels = c("V1", "V2", "V3"))
plot_data$Driver_Origin <- factor(plot_data$Driver_Origin, levels = c("V1", "V2", "V3"))
plot_data$Status <- factor(plot_data$Status, levels = c("healed", "not_healed"))

# (Da qui in poi prosegui con la Sezione 5 Plotting come prima...)
# Definiamo i 3 colori per le firme molecolari (Color blind friendly palette)

# --- SOSTITUISCI LA PARTE DEL PLOT CON QUESTA ---

p <- ggplot(plot_data, aes(x = Visit, y = Mean, group = interaction(Driver_Origin, Status))) +
  
  # A. Barre di errore (Leggere)
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE, color = Driver_Origin), 
                width = 0.1, alpha = 0.3) +
  
  # B. Linee: Usiamo alpha (trasparenza) e size (spessore) per distinguere lo Status
  geom_line(aes(color = Driver_Origin, alpha = Status, linewidth = Status)) +
  
  # C. Punti
  geom_point(aes(color = Driver_Origin, alpha = Status), size = 3) +
  
  # D. Definizione Colori (Le 3 firme)
  scale_color_manual(values = c("V1" = "red", "V2" = "blue", "V3" = "green"),
                     name = "Visit") +
  
  # E. Definizione Trasparenza (Alpha)
  # Healed = 1 (Pieno), Not Healed = 0.4 (Trasparente)
  scale_alpha_manual(values = c("healed" = 1, "not_healed" = 0.4),
                     name = "Condition") +
  
  # F. Definizione Spessore Linea
  # Healed = Spesso, Not Healed = Sottile
  scale_linewidth_manual(values = c("healed" = 1.5, "not_healed" = 0.8),
                         name = "Condition") +
  
  theme_bw() +
  labs(
    title = "Temporal analysis - unique proteins",
    y = "Average Z-Score",
    x = "Visit"
  )

print(p)
# # Salvataggio
ggsave(file.path(metadata_dir, "temporal_analysis", "temporal_analysis_unique.pdf"), plot=p, width=10, height=8)
