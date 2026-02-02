# ==============================================================================
# SETUP AMBIENTE E LIBRERIE
# ==============================================================================
# Assicurati di avere installato il pacchetto 'readxl' per leggere i file Excel
# install.packages("readxl")

library(ggplot2)
library(dplyr)
library(readxl)
library(scales) 

# ==============================================================================
# DEFINIZIONE PATH (DA COMPILARE)
# ==============================================================================
# Inserisci qui il percorso del tuo file Excel
input_xlsx_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusters.xlsx" 

# Inserisci qui la cartella dove salvare i PDF
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/barplots" 

# ==============================================================================
# CARICAMENTO E PRE-PROCESSING DATI
# ==============================================================================
if(input_xlsx_path == "" | output_dir == "") {
  stop("Error: Please define 'input_xlsx_path' and 'output_dir' before running.")
}

message("Reading Excel file...")
deg_data <- read_excel(input_xlsx_path)

# Verifica preliminare delle colonne necessarie
required_cols <- c("avg_log2FC", "p_val_adj", "cluster")
if(!all(required_cols %in% colnames(deg_data))) {
  stop("Error: The Excel file does not contain the required columns: 'avg_log2FC', 'p_val_adj', 'cluster'.")
}

message("Filtering and summarizing data...")

plot_data <- deg_data %>%
  filter(p_val_adj < 0.05) %>%
  mutate(
    Direction = ifelse(avg_log2FC > 0, "UP", "DOWN")
    # NOTA: Non convertiamo 'cluster' in factor qui per evitare l'ordinamento errato
  ) %>%
  group_by(cluster, Direction) %>%
  summarise(Gene_Count = n(), .groups = 'drop')

# --- FIX ORDINAMENTO CLUSTER (CRUCIALE) ---
# 1. Assicuriamoci che i cluster siano numerici (rimuove il problema "0, 1, 10")
plot_data$cluster <- as.numeric(as.character(plot_data$cluster))

# 2. Ordiniamo i livelli numericamente
sorted_levels <- sort(unique(plot_data$cluster))

# 3. Applichiamo l'ordine al fattore
plot_data$cluster <- factor(plot_data$cluster, levels = sorted_levels)

# Controllo: Se non ci sono geni significativi
if(nrow(plot_data) == 0) {
  stop("Warning: No genes passed the significance threshold (p_val_adj < 0.05). Check your data.")
}

# Impostiamo l'ordine della legenda (UP/DOWN)
plot_data$Direction <- factor(plot_data$Direction, levels = c("UP", "DOWN"))

# ==============================================================================
# GENERAZIONE PLOT
# ==============================================================================
message("Generating plot...")

# Definizione colori custom
custom_colors <- c("UP" = "#FFFF00", "DOWN" = "#0432FF")

p <- ggplot(plot_data, aes(x = cluster, y = Gene_Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = custom_colors) +
  # Espande l'asse Y del 10% in alto per fare spazio ai numeri
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_bw() +
  labs(
    title = "Number of Differentially Expressed Genes (DEGs) per Cluster",
    subtitle = "Significance threshold: p_adj < 0.05",
    y = "Number of Genes",
    x = "Cluster",
    fill = "Regulation"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    text = element_text(size = 12),
    # Margini ottimizzati per evitare tagli a destra
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
  )

# Aggiunge i numeri sopra le barre
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

output_file <- file.path(output_dir, "Barplot_DEGs_UP_DOWN_per_Cluster_Sorted.pdf")

# Calcolo larghezza dinamica: 
# Aumentiamo leggermente il moltiplicatore (da 0.8 a 0.9) per dare più respiro orizzontale
plot_width <- max(10, length(unique(plot_data$cluster)) * 0.9)

ggsave(filename = output_file, plot = p, width = plot_width, height = 6, dpi = 300)

message("Success! Plot saved to: ", output_file)