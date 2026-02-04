# ==============================================================================
# SETUP AMBIENTE E LIBRERIE
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales) 

# ==============================================================================
# DEFINIZIONE PATH
# ==============================================================================
input_rds_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/barplots/"

# input_rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"
# output_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/barplots/" 

# ==============================================================================
# DEFINIZIONE ANNOTAZIONE CELLULE
# ==============================================================================
cluster_names <- c(
  "0" = "T cells (cytotoxic gamma delta)", 
  "1" = "CD4 T cells (CD4 memory/activated Th2 cells)",
  "2" = "CD4 T cells (naive/CM)",
  "3" = "NK cells (CD16)",
  "4" = "CD4 T cells (activated CD4 memory T cells)",
  "5" = "CD8 T cells (EM)",
  "6" = "CD4 T cells (naive/CM)",
  "7" = "T cells (cytotoxic gamma delta)",
  "8" = "CD4 T cells (naive/CM)",
  "9" = "B cells (activated/pre-plasmablasts)", 
  "10" = "Classical monocytes",
  "11" = "B cells (naive/transitional)",
  "12" = "T cells (NKT-like gamma delta)",
  "13" = "CD4/CD8 T cells (memory T cells)",
  "14" = "B cells (naive/transitional)", 
  "15" = "Nonclassical monocytes",
  "16" = "CD4 T cells (naive/CM)",
  "17" = "Intermediate monocytes",
  "18" = "pDCs & cDC2"
)

# Formattazione con 'newline' prima della parentesi
formatted_labels <- setNames(gsub(" \\(", "\n(", cluster_names), names(cluster_names))

# ==============================================================================
# CARICAMENTO DATI E PREPARAZIONE
# ==============================================================================
if(input_rds_path == "" | output_dir == "") {
  stop("Errore: Definisci 'input_rds_path' e 'output_dir' prima di eseguire.")
}

message("Caricamento oggetto Seurat in corso...")
so <- readRDS(input_rds_path)

required_cols <- c("RNA_snn_res.0.7", "condition")
if(!all(required_cols %in% colnames(so@meta.data))) {
  stop("Errore: Le colonne richieste non esistono nei metadata.")
}

message("Creazione dataframe riassuntivo...")
df_plot <- so@meta.data %>%
  select(Cluster = RNA_snn_res.0.7, Condition = condition) %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n(), .groups = 'drop')

# Assicuriamo che 'Cluster' sia un fattore con livelli ordinati numericamente
numeric_levels <- sort(unique(as.numeric(as.character(df_plot$Cluster))))
df_plot$Cluster <- factor(df_plot$Cluster, levels = numeric_levels)

# ==============================================================================
# PREPARAZIONE ETICHETTE ASSE X
# ==============================================================================
cluster_totals <- df_plot %>%
  group_by(Cluster) %>%
  summarise(Total = sum(Count))

# Creiamo le etichette complesse: "Nome Cellula \n(n=...)"
# Usiamo formatted_labels[as.character(x)] per mappare l'ID al nome
complex_x_labels <- setNames(
  paste0(formatted_labels[as.character(cluster_totals$Cluster)]), 
  cluster_totals$Cluster
)

# ==============================================================================
# TEMA COMUNE (No Grid)
# ==============================================================================
common_theme <- theme_bw() +
  theme(
    # RIMOZIONE COMPLETA GRIGLIA
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
    legend.position = "top",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")
  )

# ==============================================================================
# PLOTS
# ==============================================================================

# --- VARIANTE 1: STACKED ---
p1 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("not_healed" = "#fdae61", "healed" = "#7fbc41")) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_x_discrete(labels = complex_x_labels) + # Etichette con conteggi
  labs(
    title = "Distribution of cells per Cell Type",
    y = "Number of cells",
    x = "Cell Type (Total count)",
    fill = "Condition"
  ) +
  common_theme

# --- VARIANTE 2: SPLIT (Faceted) ---
p2 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) + 
  facet_wrap(~Condition, scales = "fixed", ncol = 1) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_x_discrete(labels = formatted_labels) + # Etichette semplici formattate
  labs(
    title = "Cell counts per condition",
    subtitle = "Comparison: Healed vs Not Healed",
    y = "Number of cells",
    x = "Cell Type",
    fill = "Cell Type"
  ) +
  common_theme +
  theme(legend.position = "none") # Legenda ridondante qui

# --- VARIANTE 3: SIDE-BY-SIDE ---
p3 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("not_healed" = "#fdae61", "healed" = "#7fbc41")) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_x_discrete(labels = formatted_labels) + # Etichette semplici formattate
  labs(
    title = "Cell counts by Cell Type and condition",
    subtitle = "Side-by-side Comparison",
    y = "Number of cells",
    x = "Cell Type",
    fill = "Condition"
  ) +
  common_theme

# ==============================================================================
# SALVATAGGIO OUTPUT
# ==============================================================================
message("Saving plots...")

# Nomi file aggiornati "Cells"
file_stacked <- file.path(output_dir, "Barplot_Cells_Stacked_Condition.pdf")
file_split   <- file.path(output_dir, "Barplot_Cells_Split_Condition.pdf")
file_dodged  <- file.path(output_dir, "Barplot_Cells_SideBySide_Condition.pdf")

# Aumento width/height per accomodare le etichette lunghe ruotate
ggsave(filename = file_stacked, plot = p1, width = 14, height = 9, dpi = 300)
ggsave(filename = file_split, plot = p2, width = 12, height = 12, dpi = 300)
ggsave(filename = file_dodged, plot = p3, width = 14, height = 9, dpi = 300)

message("Done! Plots have been saved to: ", output_dir)