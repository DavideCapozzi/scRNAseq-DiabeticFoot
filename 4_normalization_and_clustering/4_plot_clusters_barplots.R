# ==============================================================================
# SETUP AMBIENTE E LIBRERIE
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales) # Per formattare i numeri degli assi (es. 10,000 invece di 10000)

# ==============================================================================
# DEFINIZIONE PATH (DA COMPILARE)
# ==============================================================================
# Inserisci qui il percorso del tuo file RDS
input_rds_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/seurat_res_0.7/seurat_res_0.7.rds"

# Inserisci qui la cartella dove salvare i PDF
output_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/barplots/" 

# ==============================================================================
# CARICAMENTO DATI E PREPARAZIONE
# ==============================================================================
# Verifica che i path non siano vuoti
if(input_rds_path == "" | output_dir == "") {
  stop("Errore: Definisci 'input_rds_path' e 'output_dir' prima di eseguire.")
}

message("Caricamento oggetto Seurat in corso...")
so <- readRDS(input_rds_path)

# Verifica esistenza colonne
required_cols <- c("RNA_snn_res.0.7", "condition")
if(!all(required_cols %in% colnames(so@meta.data))) {
  stop("Errore: Le colonne 'RNA_snn_res.0.7' o 'condition' non esistono nei metadata.")
}

message("Creazione dataframe riassuntivo...")
# Estrazione metadata e creazione tabella di frequenza
# Best practice: Aggregare prima di plottare per gestire grandi numeri di cellule
df_plot <- so@meta.data %>%
  select(Cluster = RNA_snn_res.0.7, Condition = condition) %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n(), .groups = 'drop')

# Assicuriamoci che i cluster siano ordinati correttamente (essendo Factor dovrebbe essere ok)
# Se necessario, forza l'ordine dei livelli qui.

# ==============================================================================
# PREPARAZIONE ETICHETTE ASSE X (CONTE TOTALI)
# ==============================================================================
# Calcoliamo i totali per cluster per visualizzarli nell'asse X del primo grafico
cluster_totals <- df_plot %>%
  group_by(Cluster) %>%
  summarise(Total = sum(Count))

# Creiamo un vettore nominato per mappare: "ID Cluster" -> "ID Cluster \n (n=...)"
# Questo ci permette di modificare le etichette dell'asse X senza toccare i dati originali
x_axis_labels <- setNames(
  paste0(cluster_totals$Cluster, "\n(n=", comma(cluster_totals$Total), ")"), 
  cluster_totals$Cluster
)

# ==============================================================================
# VARIANTE 1: STACKED BARPLOT (Una sopra l'altra)
# ==============================================================================
# Shows the absolute proportion of conditions for each cluster with total counts on X

p1 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("not_healed" = "#E69F00", "healed" = "#56B4E9")) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  # Qui applichiamo le etichette personalizzate con i numeri totali
  scale_x_discrete(labels = x_axis_labels) + 
  theme_bw() +
  labs(
    title = "Distribution of cells per cluster",
    subtitle = "Clustering Resolution: 0.7",
    y = "Number of cells",
    x = "Cluster (Total count)",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9), # Ruota e ridimensiona per leggibilità
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# ==============================================================================
# VARIANTE 2: SPLIT BARPLOT (Separati per Condizione)
# ==============================================================================
# Comparison of distribution shapes between conditions (Faceted)

p2 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) + # 'linewidth' sostituisce 'size' nelle nuove versioni di ggplot2
  facet_wrap(~Condition, scales = "fixed", ncol = 1) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  labs(
    title = "Cell counts per condition",
    subtitle = "Comparison: Healed vs Not Healed",
    y = "Number of cells",
    x = "Cluster",
    fill = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "none" 
  )

# ==============================================================================
# VARIANTE 3: SIDE-BY-SIDE BARPLOT (Una accanto all'altra)
# ==============================================================================
# Direct comparison of counts per condition within each cluster

p3 <- ggplot(df_plot, aes(x = Cluster, y = Count, fill = Condition)) +
  # position_dodge mette le barre una accanto all'altra
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("not_healed" = "#E69F00", "healed" = "#56B4E9")) + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  labs(
    title = "Cell counts by cluster and condition",
    subtitle = "Side-by-side Comparison",
    y = "Number of cells",
    x = "Cluster",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# ==============================================================================
# SALVATAGGIO OUTPUT
# ==============================================================================
message("Saving plots...")

# Definisci nomi file
file_stacked <- file.path(output_dir, "Barplot_Clusters_Stacked_Condition.pdf")
file_split   <- file.path(output_dir, "Barplot_Clusters_Split_Condition.pdf")
file_dodged  <- file.path(output_dir, "Barplot_Clusters_SideBySide_Condition.pdf")

# Salva PDF 1 (Stacked con conteggi su asse X)
ggsave(filename = file_stacked, plot = p1, width = 12, height = 7, dpi = 300)

# Salva PDF 2 (Split verticalmente)
ggsave(filename = file_split, plot = p2, width = 10, height = 9, dpi = 300)

# Salva PDF 3 (Side-by-side)
ggsave(filename = file_dodged, plot = p3, width = 12, height = 7, dpi = 300)

message("Done! Plots have been saved to: ", output_dir)