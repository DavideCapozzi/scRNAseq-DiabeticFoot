# Installa e carica il pacchetto richiesto (se non l'hai già fatto)
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

# --- Impostazioni ---

# Directory base fornita (la parte che contiene le sottocartelle V1, V2, V3)
base_dir <- "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"

# Nomi delle visite/punti temporali
visite <- c("V1", "V2", "V3")

# Geni di interesse
geni_target <- c("CCL4", "CHL1", "COMP")

# --- Caricamento e Preparazione dei Dati ---

# Inizializza un data frame vuoto per raccogliere i dati
dati_combinati <- data.frame()

# Cicla attraverso ogni visita per caricare il file DEG e filtrare i geni
for (visita in visite) {
  # Costruisci il percorso del file DEG.txt per la visita corrente
  percorso_file <- file.path(base_dir, paste0("DiabeticFoot_", visita), paste0("Results_", visita), "DEGs", "DEG.txt")
  
  # Controlla se il file esiste prima di provare a leggerlo
  if (file.exists(percorso_file)) {
    
    # Leggi il file DEG (presuppone che sia delimitato da tab)
    # Si noti che il file allegato usa 'gene' come nome della prima colonna
    deg_data <- read.delim(percorso_file, header = TRUE, sep = "\t")
    
    # Filtra per i geni target
    dati_filtrati <- deg_data %>%
      filter(gene %in% geni_target) %>%
      select(gene, log2.FC) # Seleziona solo il gene e il log2-FC
    
    # Aggiungi la colonna della visita (punto temporale)
    dati_filtrati$visita <- visita
    
    # Aggiungi i dati filtrati al data frame combinato
    dati_combinati <- bind_rows(dati_combinati, dati_filtrati)
    
  } else {
    warning(paste("File non trovato per la visita:", visita, "al percorso:", percorso_file))
  }
}

# Assicurati che 'visita' sia un fattore ordinato per un corretto ordinamento sull'asse X
dati_combinati$visita_ordinata <- factor(dati_combinati$visita, levels = visite)

# Rinomina la colonna log2.FC per chiarezza nel grafico
dati_combinati <- dati_combinati %>%
  rename(LogFC = log2.FC)

# --- Generazione del Grafico ---
# 
grafico_logfc <- ggplot(dati_combinati, aes(x = visita_ordinata, y = LogFC, group = gene, color = gene)) +
  # Aggiungi linee e punti
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  
  # Personalizzazione dei titoli e delle etichette
  labs(
    title = "LogFC (Healed vs Not healed)",
    x = "Visit",
    y = expression(LogFC),
    color = "Signature Protein" # Titolo della legenda per il colore
  ) +
  
  # Personalizzazione del tema per un aspetto pulito (opzionale)
  theme_minimal() +
  
  # Personalizzazione dei colori (se necessario, altrimenti usa i colori di default di ggplot)
  scale_color_manual(values = c("CCL4" = "#E41A1C", "CHL1" = "#377EB8", "COMP" = "#4DAF4A")) +
  
  # Regola l'asse Y per replicare l'intervallo del grafico di esempio (circa da -1 a 2)
  # Usa coord_cartesian per zoomare senza rimuovere i dati
  coord_cartesian(ylim = c(min(dati_combinati$LogFC) - 0.2, max(dati_combinati$LogFC) + 0.2)) +
  
  # Opzioni aggiuntive per il testo
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Stampa il grafico
print(grafico_logfc)

# Se vuoi salvare il grafico su un file:
ggsave("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/temporal_analysis/trajectories_logFC.png", plot = grafico_logfc, width = 8, height = 6, dpi = 300)
