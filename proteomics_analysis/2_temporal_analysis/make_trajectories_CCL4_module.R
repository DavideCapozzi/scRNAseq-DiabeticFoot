library(tidyverse)

# --- 1. Configurazione Percorsi ---
base_dir <- "/Users/federicadannunzio/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"

# File dei nodi per identificare il modulo (da V1)
path_nodes <- file.path(base_dir, "DiabeticFoot_V1", "Results_V1", "txtFile", "cytoscapeInput_nodes.txt")

# --- 2. Identificazione Geni del Modulo ---
nodes_data <- read.delim(path_nodes, header = TRUE, sep = "\t", check.names = FALSE)

# Identifichiamo il colore del modulo di CCL4 (colonna 3) [cite: 1]
target_module <- nodes_data %>%
  filter(nodeName == "CCL4") %>%
  pull(3) 

# Estraiamo tutti i geni appartenenti a questo modulo [cite: 1]
geni_modulo <- nodes_data %>%
  filter(nodes_data[[3]] == target_module) %>%
  pull(nodeName)

# --- 3. Caricamento Dati ---
visite <- c("V1", "V2", "V3")
dati_combinati <- data.frame()

for (v in visite) {
  percorso_file <- file.path(base_dir, paste0("DiabeticFoot_", v), paste0("Results_", v), "DEGs", "DEG.txt")
  
  if (file.exists(percorso_file)) {
    # check.names = FALSE per gestire "log2-FC" con il trattino [cite: 2]
    deg_data <- read.delim(percorso_file, header = TRUE, sep = "\t", check.names = FALSE)
    
    dati_filtrati <- deg_data %>%
      filter(gene %in% geni_modulo) %>%
      mutate(
        visita = v,
        # Etichetta asterisco se p-value < 0.05 [cite: 2]
        label_sig = ifelse(pval < 0.05, "*", "")
      ) %>%
      select(gene, `log2-FC`, visita, pval, label_sig)
    
    dati_combinati <- bind_rows(dati_combinati, dati_filtrati)
  }
}

dati_combinati$visita <- factor(dati_combinati$visita, levels = visite)

# --- 4. Plot delle Traiettorie ---
grafico_final <- ggplot(dati_combinati, aes(x = visita, y = `log2-FC`, group = gene, color = gene)) +
  # Linee spesse per massima visibilità
  geom_line(linewidth = 1.5, alpha = 0.8) +
  
  # Pallini leggermente più grandi per contenere l'asterisco
  geom_point(size = 5) + 
  
  # Asterisco sovrapposto esattamente al pallino
  geom_text(
    aes(label = label_sig), 
    color = "black",      # Bianco per contrasto dentro il pallino colorato
    size = 10,            # Dimensione proporzionata al pallino
    fontface = "bold",
    hjust = 0.5,          # Centratura orizzontale perfetta
    vjust = 0.6,         # Centratura verticale (corregge l'offset naturale del carattere *)
    show.legend = FALSE
  ) +
  
  labs(
    title = "LogFC trajectories - CCL4 module",
    x = "Visit",
    y = expression(Log[2]*" FC"),
    color = "Proteins"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid.major = element_line(color = "grey90")
  )

print(grafico_final)

