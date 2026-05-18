library(tidyverse)
library(xlsx)

base_dir <- "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset"
visite <- c("V1", "V2", "V3")
geni_target <- c("CCL4", "CHL1", "COMP")

dati_combinati <- data.frame()

for (visita in visite) {
  percorso_file <- file.path(base_dir, paste0("DiabeticFoot_", visita), paste0("Results_", visita), "DEGs", "DEG.txt")
  
  if (file.exists(percorso_file)) {
    deg_data <- read.delim(percorso_file, header = TRUE, sep = "\t")
    
    dati_filtrati <- deg_data %>%
      filter(gene %in% geni_target) %>%
      select(gene, log2.FC, pval)
    
    dati_filtrati$visita <- visita
    dati_combinati <- bind_rows(dati_combinati, dati_filtrati)
  } else {
    warning(paste("File non trovato per la visita:", visita, "al percorso:", percorso_file))
  }
}

dati_combinati$visita_ordinata <- factor(dati_combinati$visita, levels = visite)
dati_combinati <- dati_combinati %>%
  rename(LogFC = log2.FC)

grafico_logfc_pannelli <- ggplot(dati_combinati, aes(x = visita_ordinata, y = LogFC, group = gene, color = gene)) +
  geom_line(linewidth = 0.5, lineend = "round",
            linejoin = "round") +
  #geom_smooth(method = "loess", se = FALSE, span = 0.7, linewidth = 0.5) +
  geom_point(size = 1.5) +
  labs(
    title = "LogFC Trajectories (Healed vs Not healed)",
    x = "Visit",
    y = expression(LogFC),
    color = "Signature Protein"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("CCL4" = "green", "CHL1" = "blue", "COMP" = "turquoise")) +
  facet_wrap(~ gene, ncol = length(geni_target), scales = "fixed") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor = element_line(color = "grey90", linetype = "dotted"),
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(grafico_logfc_pannelli)

###########################
ggsave("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/temporal_analysis/trajectories_logFC.png", plot = grafico_logfc_pannelli, width = 10, height = 4, dpi = 300)
write.xlsx(dati_combinati, "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/temporal_analysis/signature_proteins_pval_fc.xlsx")

##########################
curve_df <- dati_combinati %>%
  group_by(gene) %>%
  arrange(visita_ordinata) %>%
  summarise(
    x = list(as.numeric(visita_ordinata)),
    y = list(LogFC),
    .groups = "drop"
  ) %>%
  mutate(
    curve = map2(x, y, ~ {
      xs <- seq(min(.x), max(.x), length.out = 100)
      data.frame(
        visita_num = xs,
        LogFC = spline(.x, .y, xout = xs)$y
      )
    })
  ) %>%
  unnest(curve)

p_smooth <- p_smooth +
  theme(
    # Sfondo pannelli
    panel.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    
    # Griglia puntinata
    panel.grid.major = element_line(
      color = "grey75",
      linetype = "dotted",
      linewidth = 0.6
    ),
    panel.grid.minor = element_line(
      color = "grey85",
      linetype = "dotted",
      linewidth = 0.4
    ),
    
    # Assi
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    
    # Testi
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    # Facet strip
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    
    # Niente legenda (coerente col plot)
    legend.position = "none"
  )

p_smooth <- p_smooth +
  scale_color_manual(
    values = c(
      "CCL4" = "green",  
      "CHL1" = "blue",  
      "COMP" = "turquoise"   
    )
  )


print(p_smooth)

ggsave("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/temporal_analysis/trajectories_logFC_smooth.png", plot = p_smooth, width = 10, height = 4, dpi = 300)
