# ==============================================================================
# SCRIPT WGCNA: TRAIETTORIE UNIFICATE (Senza distinzione di Rank)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)

options(stringsAsFactors = FALSE)

# 1. PERCORSI
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

# 2. CARICAMENTO METADATI
# ==============================================================================
if (!file.exists(metadata_file)) stop("Metadati non trovati")
metadata <- read_excel(metadata_file)

meta_clean <- metadata %>%
  select(`Proteomics ID`, condition) %>%
  rename(SampleID_Clean = `Proteomics ID`, Status = condition) %>%
  mutate(SampleID_Clean = trimws(SampleID_Clean))

# 3. FUNZIONI HELPER
# ==============================================================================

get_top_hub_genes <- function(result_path, visit_name) {
  p1 <- file.path(result_path, "txtFile", "geneInfo.txt")
  p2 <- file.path(result_path, "geneInfo.txt")
  if (file.exists(p1)) f <- p1 else f <- p2
  
  geneInfo <- read.table(f, header = TRUE, sep = "\t")
  
  if(colnames(geneInfo)[1] == "moduleColor") {
    geneInfo <- geneInfo %>% rownames_to_column(var = "GeneName")
  } else {
    colnames(geneInfo)[1] <- "GeneName"
  }
  
  module_stats <- geneInfo %>%
    filter(moduleColor != "grey") %>%
    group_by(moduleColor) %>%
    summarise(mean_GS = mean(abs(get(paste0("GS.", trait_col_wgcna))), na.rm=T)) %>%
    arrange(desc(mean_GS))
  
  if(nrow(module_stats)==0) return(NULL)
  top_mod <- module_stats$moduleColor[1]
  
  mm_col <- paste0("MM.", top_mod)
  gs_col <- paste0("GS.", trait_col_wgcna)
  
  top_genes_df <- geneInfo %>%
    filter(moduleColor == top_mod) %>%
    mutate(HubScore = abs(get(gs_col)) + abs(get(mm_col))) %>%
    arrange(desc(HubScore)) %>%
    distinct(GeneName, .keep_all = TRUE) %>%
    head(5) %>% 
    select(GeneName)
  
  top_genes_df$Driver_Origin <- visit_name
  top_genes_df$Module <- top_mod
  
  return(top_genes_df)
}

get_expression_data <- function(result_path, visit_name, target_genes, metadata_df) {
  r1 <- file.path(result_path, "RData", "dataInput.RData")
  r2 <- file.path(result_path, "dataInput.RData")
  if (file.exists(r1)) rf <- r1 else rf <- r2
  
  env <- new.env()
  load(rf, envir = env)
  datExpr <- env$datExpr
  
  dataset_cols <- colnames(datExpr)
  valid_genes <- c()
  
  for(g in target_genes) {
    if(g %in% dataset_cols) {
      valid_genes <- c(valid_genes, g)
    } else {
      g_alt <- gsub("-", ".", g)
      if(g_alt %in% dataset_cols) {
        colnames(datExpr)[which(colnames(datExpr) == g_alt)] <- g
        valid_genes <- c(valid_genes, g)
      } else {
        g_alt2 <- gsub("\\.", "-", g)
        if(g_alt2 %in% dataset_cols) {
          colnames(datExpr)[which(colnames(datExpr) == g_alt2)] <- g
          valid_genes <- c(valid_genes, g)
        }
      }
    }
  }
  
  if(length(valid_genes) == 0) return(NULL)
  
  expr_subset <- datExpr[, valid_genes, drop = FALSE]
  
  df_long <- as.data.frame(expr_subset) %>%
    rownames_to_column("SampleID_Original") %>%
    pivot_longer(cols = all_of(valid_genes), names_to = "GeneName", values_to = "Expression") %>%
    mutate(Visit = visit_name)
  
  df_long$Match_Key <- tolower(gsub("\\s+", "", df_long$SampleID_Original))
  metadata_df$Match_Key <- tolower(gsub("\\s+", "", metadata_df$SampleID_Clean))
  
  df_merged <- df_long %>%
    left_join(metadata_df, by = "Match_Key") %>%
    filter(!is.na(Status))
  
  return(df_merged)
}

# 4. ESECUZIONE
# ==============================================================================

# A. Trova Geni
all_targets <- data.frame()
print("Ricerca geni...")
for(v in visits) {
  tmp <- get_top_hub_genes(visit_paths[[v]], v)
  if(!is.null(tmp)) all_targets <- rbind(all_targets, tmp)
}
unique_genes <- unique(all_targets$GeneName)
print(paste("Geni trovati:", length(unique_genes)))

# B. Estrai Dati
full_data <- data.frame()
print("Estrazione dati...")
for(v in visits) {
  tmp <- get_expression_data(visit_paths[[v]], v, unique_genes, meta_clean)
  if(!is.null(tmp)) full_data <- rbind(full_data, tmp)
}

if(nrow(full_data) == 0) stop("Errore: Dataset finale vuoto.")

# C. Calcolo Z-Score e Statistiche
plot_data_prep <- full_data %>%
  inner_join(all_targets, by = "GeneName") 

plot_data <- plot_data_prep %>%
  group_by(GeneName) %>%
  mutate(Z_Score = scale(Expression)) %>%
  ungroup() %>%
  group_by(GeneName, Visit, Status, Driver_Origin) %>%
  summarise(
    Mean_Z = mean(Z_Score, na.rm=T), 
    SE_Z = sd(Z_Score, na.rm=T)/sqrt(n()), 
    .groups="drop"
  )

plot_data$Visit <- factor(plot_data$Visit, levels=c("V1","V2","V3"))

# Etichette Pannelli
plot_data$Panel_Label <- paste("Drivers identified in", plot_data$Driver_Origin)

# 5. PLOT PULITO (Senza Rank)
# ==============================================================================

p <- ggplot(plot_data, aes(x = Visit, y = Mean_Z, group = interaction(GeneName, Status))) +
  
  geom_hline(yintercept = 0, linetype="solid", color="grey85", linewidth=0.5) +
  
  # Linee: Solo colore per status, tutte solide
  geom_line(aes(color = Status), linewidth = 1, alpha = 0.8) +
  
  # Punti: Solo colore per status
  geom_point(aes(color = Status), size = 2.5) +
  
  # Error bars
  geom_errorbar(aes(ymin = Mean_Z - SE_Z, ymax = Mean_Z + SE_Z, color = Status), width = 0.1, alpha = 0.6) +
  
  # Pannelli per visita di origine
  facet_wrap(~Panel_Label, ncol = 3) +
  
  # Colori
  scale_color_manual(values = c("healed" = "#009E73", "not_healed" = "#D55E00"), 
                     labels = c("Healed", "Not Healed")) +
  
  theme_bw() +
  
  labs(
    title = "Evoluzione Temporale dei Driver Molecolari (Per Visita)", 
    subtitle = "Z-Score medio dei top 3 geni identificati in ciascuna visita.",
    y = "Z-Score Expression (Mean +/- SE)", 
    x = "Time Point",
    color = "Clinical Outcome"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank() # Pulisce il background
  )

print(p)

# output_file <- file.path(metadata_dir, "Trajectory_Clean_NoRank.pdf")
# ggsave(output_file, plot = p, width = 14, height = 8)
# print(paste("Grafico salvato in:", output_file))