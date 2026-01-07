# LIBRERIE
# ==========================
library(ggplot2)
library(dplyr)

# ==========================
# DATI
# ==========================
contamination_data <- data.frame(
  sample = c("Patient1", "Patient2", "Patient3", "Patient4",
             "Patient5", "Patient6", "Patient7", "Patient8",
             "Patient9", "Patient10", "Patient11", "Patient12"),
  contamination_pct = c(2.5, 2.1, 4.5, 4.9, 9.5, 3.5, 2.6, 1.6, 1.5, 2.2, 4.7, 1.3)
)

# Converti la colonna sample in factor con l'ordine corretto
contamination_data$sample <- factor(contamination_data$sample, 
                                    levels = c("Patient1", "Patient2", "Patient3", "Patient4",
                                               "Patient5", "Patient6", "Patient7", "Patient8",
                                               "Patient9", "Patient10", "Patient11", "Patient12"))

# ==========================
# GRAFICO
# ==========================
# Barplot della contaminazione per campione
ggplot(contamination_data, aes(x = sample, y = contamination_pct, fill = sample)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(contamination_pct, "%")), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)  # Aggiunge gli assi neri
  ) +
  labs(
    title = "Contamination percentage per sample",
    x = "Samples",
    y = "Contamination (%)",
    fill = "Sample"
  ) +
  ylim(0, max(contamination_data$contamination_pct) + 2)

# ==========================
# SALVATAGGIO
# ==========================
ggsave("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/01_soupx_analysis/contamination_per_sample.png", width = 10, height = 6)