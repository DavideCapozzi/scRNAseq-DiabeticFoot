library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# =============================================================================
# CONFIGURAZIONE PERSONALIZZABILE
# =============================================================================

# Directory di output
output_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/2_pre_filtering/doublet_barplots/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir) 
} else {
  print("the dir already exists")
}

# Colori personalizzabili
colors_condition <- c(
  "healed" = "#7fbc41",
  "not_healed" = "#f46d43"
)

colors_stage <- c(
  "Before filtering" = "#252525",
  "After filtering" = "#d9d9d9"
)

# Dimensioni grafici
plot_width <- 15
plot_height <- 10
dpi <- 300

# =============================================================================
# DATI DAI LOG - ESTRATTI DIRETTAMENTE DAL remove_doublets_708.log
# =============================================================================

doublet_data <- data.frame(
  fastq_file_id = c(
    "P33554_1001", "P33554_1002", "P33554_1003", "P33554_1004",
    "P34304_1001", "P34304_1002", "P34304_1003", "P34304_1004",
    "P34304_1005", "P34304_1007", "P34304_1008", "P34304_1009"
  ),
  initial_cells = c(
    40551, 22248, 34591, 38555,
    53988, 25163, 28622, 29434,
    22305, 26392, 26323, 21452
  ),
  doublets = c(
    5714, 2252, 3704, 4610,
    6299, 2474, 3227, 3536,
    1940, 2596, 1820, 2090
  ),
  doublet_rate = c(
    14.09, 10.12, 10.71, 11.96,
    11.67, 9.83, 11.27, 12.01,
    8.70, 9.84, 6.91, 9.74
  ),
  outcome = c(
    "not_healed", "not_healed", "healed", "healed",
    "not_healed", "healed", "healed", "healed",
    "not_healed", "not_healed", "healed", "healed"
  ),
  stringsAsFactors = FALSE
)

# Tabella di mapping FastQ ID -> Patient Name
mapping_table <- data.frame(
  fastq_file_id = c(
    "P33554_1001", "P33554_1002", "P33554_1003", "P33554_1004",
    "P34304_1001", "P34304_1002", "P34304_1003", "P34304_1004",
    "P34304_1005", "P34304_1007", "P34304_1008", "P34304_1009"
  ),
  patient_name = c(
    "Patient1", "Patient2", "Patient3", "Patient4",
    "Patient5", "Patient6", "Patient7", "Patient8",
    "Patient9", "Patient10", "Patient11", "Patient12"
  ),
  stringsAsFactors = FALSE
)

# Join con la tabella di mapping
doublet_data <- doublet_data %>%
  left_join(mapping_table, by = "fastq_file_id")

# Ordinare per Patient (numericamente)
doublet_data <- doublet_data %>%
  mutate(
    patient_num = as.numeric(sub("Patient", "", patient_name))
  ) %>%
  arrange(patient_num) %>%
  select(-patient_num)

# Calcola cellule finali e aggiungi condizione
doublet_data <- doublet_data %>%
  mutate(
    final_cells = initial_cells - doublets,
    condition_label = ifelse(outcome == "healed", "Healed", "Not Healed")
  )

# Verifica totali
cat("=== SUMMARY ===\n")
cat("Total samples:", nrow(doublet_data), "\n")
cat("Total initial cells:", sum(doublet_data$initial_cells), "\n")
cat("Total doublets removed:", sum(doublet_data$doublets), "\n")
cat("Total final cells:", sum(doublet_data$final_cells), "\n")
cat("Overall doublet rate:", 
    round(sum(doublet_data$doublets) / sum(doublet_data$initial_cells) * 100, 2), "%\n\n")

# =============================================================================
# GRAFICO 1: BARPLOT AFFIANCATO (Before/After per ogni campione)
# =============================================================================

plot_data <- doublet_data %>%
  select(patient_name, condition_label, initial_cells, final_cells) %>%
  pivot_longer(
    cols = c(initial_cells, final_cells),
    names_to = "stage",
    values_to = "cells"
  ) %>%
  mutate(
    stage = factor(stage, 
                   levels = c("initial_cells", "final_cells"),
                   labels = c("Before filtering", "After filtering"))
  )

p1 <- ggplot(plot_data, 
             aes(x = factor(patient_name, levels = doublet_data$patient_name), 
                 y = cells, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = comma(cells)), 
            position = position_dodge(width = 0.7),
            vjust = -0.5,
            size = 3) +
  scale_fill_manual(values = colors_stage) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Cell count before and after doublets removal",
    subtitle = paste0("Total: ", format(sum(doublet_data$initial_cells), big.mark = ","), 
                      " → ", format(sum(doublet_data$final_cells), big.mark = ","), " cells"),
    x = "Patient",
    y = "Number of cells",
    fill = "Stage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "top"
  )

print(p1)

ggsave(paste0(output_dir, "01_cells_before_after_by_patient.pdf"), 
       plot = p1, width = plot_width, height = plot_height, dpi = dpi)

# =============================================================================
# SALVA STATISTICHE
# =============================================================================

stats_summary <- doublet_data %>%
  group_by(condition_label) %>%
  summarise(
    n_samples = n(),
    total_initial_cells = sum(initial_cells),
    total_final_cells = sum(final_cells),
    total_doublets = sum(doublets),
    mean_doublet_rate = round(mean(doublet_rate), 2),
    median_doublet_rate = round(median(doublet_rate), 2),
    min_doublet_rate = round(min(doublet_rate), 2),
    max_doublet_rate = round(max(doublet_rate), 2),
    .groups = "drop"
  )

write.csv(doublet_data, 
          paste0(output_dir, "doublet_filtering_data.csv"), 
          row.names = FALSE)
write.csv(stats_summary, 
          paste0(output_dir, "summary_statistics.csv"), 
          row.names = FALSE)

cat("\n=== STATISTICHE PER OUTCOME ===\n")
print(stats_summary)
cat("\n=== DATI DETTAGLIATI ===\n")
print(doublet_data)

