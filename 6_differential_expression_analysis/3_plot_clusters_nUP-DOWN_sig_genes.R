# ==============================================================================
# SETUP AMBIENTE E LIBRERIE
# ==============================================================================
library(ggplot2)
library(dplyr)
library(readxl)
library(scales) 

# ==============================================================================
# DEFINIZIONE PATH
# ==============================================================================
input_xlsx_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusters.xlsx" 
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/barplots" 

# ==============================================================================
# DEFINIZIONE ANNOTAZIONE CELLULE
# ==============================================================================
# Mappatura ID Cluster -> Nome Cellula
# L'ordine è preservato dall'indice (0, 1, 2...)
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

# Formattazione per andare a capo prima della parentesi
# Sostituisce " (" con "\n(" per migliorare la leggibilità nel plot
formatted_labels <- setNames(gsub(" \\(", "\n(", cluster_names), names(cluster_names))

# ==============================================================================
# CARICAMENTO E PRE-PROCESSING DATI
# ==============================================================================
if(input_xlsx_path == "" | output_dir == "") {
  stop("Error: Please define 'input_xlsx_path' and 'output_dir' before running.")
}

message("Reading Excel file...")
deg_data <- read_excel(input_xlsx_path)

required_cols <- c("avg_log2FC", "p_val_adj", "cluster")
if(!all(required_cols %in% colnames(deg_data))) {
  stop("Error: The Excel file does not contain the required columns.")
}

message("Filtering and summarizing data...")

plot_data <- deg_data %>%
  filter(p_val_adj < 0.05) %>%
  mutate(
    Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")
  ) %>%
  group_by(cluster, Direction) %>%
  summarise(Gene_Count = n(), .groups = 'drop')

# --- GESTIONE ORDINAMENTO CLUSTER ---
# Convertiamo in numerico per ordinare correttamente (0, 1, 2... 10)
plot_data$cluster <- as.numeric(as.character(plot_data$cluster))
sorted_levels <- sort(unique(plot_data$cluster))
plot_data$cluster <- factor(plot_data$cluster, levels = sorted_levels)

if(nrow(plot_data) == 0) {
  stop("Warning: No genes passed the significance threshold.")
}

plot_data$Direction <- factor(plot_data$Direction, levels = c("UP", "DOWN"))

# ==============================================================================
# GENERAZIONE PLOT
# ==============================================================================
message("Generating plot...")

custom_colors <- c("UP" = "#FFFF00", "DOWN" = "#0432FF")

p <- ggplot(plot_data, aes(x = cluster, y = Gene_Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  # Qui applichiamo la mappatura delle etichette mantenendo l'ordine dei livelli numerici
  scale_x_discrete(labels = formatted_labels) +
  theme_bw() +
  labs(
    title = "Number of Differentially Expressed Genes (DEGs) per Cell Type",
    subtitle = "Significance threshold: p_adj < 0.05",
    y = "Number of Genes",
    x = "Cell Type", # Etichetta asse aggiornata
    fill = "Regulation"
  ) +
  theme(
    # Rotazione etichette e allineamento per evitare sovrapposizioni
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), 
    # RIMOZIONE GRIGLIE (Sfondo pulito)
    panel.grid = element_blank(),
    legend.position = "top",
    text = element_text(size = 12),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 60, unit = "pt") # Margine sinistro aumentato per i nomi lunghi
  )

p <- p + geom_text(
  aes(label = Gene_Count), 
  position = position_dodge(width = 0.8), 
  vjust = -0.5, 
  size = 3
)

# ==============================================================================
# SALVATAGGIO OUTPUT
# ==============================================================================
message("Saving PDF...")

# Nome file aggiornato da Clusters a Cells
output_file <- file.path(output_dir, "Barplot_DEGs_UP_DOWN_per_Cells_Sorted.pdf")

# Calcolo larghezza dinamica aumentato per accomodare i nomi lunghi delle cellule
plot_width <- max(12, length(unique(plot_data$cluster)) * 1.2) 

ggsave(filename = output_file, plot = p, width = plot_width, height = 8, dpi = 300) # Altezza aumentata per le etichette ruotate

message("Success! Plot saved to: ", output_file)