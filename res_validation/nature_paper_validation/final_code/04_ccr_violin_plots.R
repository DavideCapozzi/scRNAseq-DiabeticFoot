# ==============================================================================
# SCRIPT 04: CCR1 / CCR5 / CCR8 — VIOLIN PLOTS
# Dataset: GSE165816 - DFU scRNA-seq
#
# SECTIONS:
#   A. Gene presence check
#   B. Violin plots — all three conditions (per cell type)
#   C. Violin plots — Healing vs Non-Healing (per cell type)
#   D. Violin plots — Diabetes subset (per cell type)
#   E. Combined panel (one PDF per receptor, tutti e tre i plot in colonna)
#
# PREREQUISITE: Script 02 completed (seurat_obj_annotated.rds present)
# ==============================================================================

rm(list = ls())
set.seed(42)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ==============================================================================
# 0. SETUP — percorsi e tema (identici allo script CCL4)
# ==============================================================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/nature_paper_validation"
res_dir  <- file.path(base_dir, "res")
plot_dir <- file.path(res_dir, "plots", "CCR")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- tema globale (identico a theme_dfu dello script CCL4) ----------
theme_dfu <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line         = element_line(color = "black", linewidth = 0.5),
      axis.ticks        = element_line(color = "black"),
      axis.text         = element_text(color = "black"),
      axis.title        = element_text(color = "black"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      strip.background  = element_rect(fill = "grey92", color = "black", linewidth = 0.4),
      strip.text        = element_text(color = "black", face = "bold"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.text       = element_text(color = "black"),
      legend.title      = element_text(color = "black", face = "bold"),
      plot.title        = element_text(face = "bold", color = "black", size = base_size + 1),
      plot.subtitle     = element_text(color = "grey35", size = base_size - 1),
      plot.caption      = element_text(color = "grey50",  size = base_size - 2)
    )
}

seurat_theme <- function() {
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black"),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    strip.background  = element_rect(fill = "grey92", color = "black"),
    strip.text        = element_text(color = "black", face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black", face = "bold"),
    plot.title        = element_text(face = "bold", color = "black"),
    plot.subtitle     = element_text(color = "grey35")
  )
}

cond_colors <- c("Healing" = "#2196F3", "Non-Healing" = "#F44336", "Diabetes" = "#4CAF50")

# ==============================================================================
# 1. CARICA OGGETTO ANNOTATO
# ==============================================================================
cat("--- Loading annotated object ---\n")
seurat_obj <- readRDS(file.path(res_dir, "seurat_obj_annotated.rds"))

cat(sprintf("Total cells: %d | Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))
cat("Cells per condition:\n"); print(table(seurat_obj$Condition))

# Palette colori per cell type (stessa logica dello script CCL4)
n_ct      <- length(unique(seurat_obj$CellType))
ct_colors <- scales::hue_pal()(n_ct)
names(ct_colors) <- sort(unique(seurat_obj$CellType))

# ==============================================================================
# SECTION A: CHECK PRESENZA DEI GENI
# ==============================================================================
receptors_target <- c("CCR1", "CCR5", "CCR8")

receptors_ok <- receptors_target[receptors_target %in% rownames(seurat_obj)]
receptors_missing <- setdiff(receptors_target, receptors_ok)

if (length(receptors_missing) > 0) {
  cat(sprintf("[WARNING] Gene/s non trovati nell'oggetto: %s\n",
              paste(receptors_missing, collapse = ", ")))
  cat("  Controlla varianti di nome con: grep('CCR', rownames(seurat_obj), value = TRUE)\n")
}

if (length(receptors_ok) == 0) {
  stop("[ERROR] Nessuno dei recettori target trovato nell'oggetto Seurat.")
}

cat(sprintf("\nRecettori da analizzare: %s\n", paste(receptors_ok, collapse = ", ")))

# Statistiche di base per ciascun recettore
cat("\n=== Baseline Statistics ===\n")
for (gene in receptors_ok) {
  expr_vec <- FetchData(seurat_obj, vars = gene)[[gene]]
  cat(sprintf(
    "  %s — positive cells: %d / %d (%.1f%%) | mean expr (all): %.4f\n",
    gene,
    sum(expr_vec > 0), length(expr_vec),
    100 * mean(expr_vec > 0),
    mean(expr_vec)
  ))
}

# ==============================================================================
# Subset per le sezioni successive
# ==============================================================================
Idents(seurat_obj) <- "CellType"

# Healing + Non-Healing (senza Diabetes)
seurat_hnh <- subset(seurat_obj, subset = Condition %in% c("Healing", "Non-Healing"))
Idents(seurat_hnh) <- "CellType"
cat(sprintf("\nHealing: %d cells | Non-Healing: %d cells\n",
            sum(seurat_hnh$Condition == "Healing"),
            sum(seurat_hnh$Condition == "Non-Healing")))

# Solo Diabetes (controllo)
seurat_diab <- subset(seurat_obj, subset = Condition == "Diabetes")
Idents(seurat_diab) <- "CellType"
cat(sprintf("Diabetes: %d cells | Samples: %s\n",
            ncol(seurat_diab),
            paste(unique(seurat_diab$orig.ident), collapse = ", ")))

# ==============================================================================
# SECTIONS B / C / D — loop per ogni recettore
# ==============================================================================

for (gene in receptors_ok) {

  cat(sprintf("\n\n========== Processing: %s ==========\n", gene))

  # --------------------------------------------------------------------------
  # SECTION B: Violin plot — tutte e tre le condizioni, per cell type
  # --------------------------------------------------------------------------
  cat(sprintf("  [B] VlnPlot all conditions — %s\n", gene))

  p_vln_all <- VlnPlot(
    seurat_obj,
    features = gene,
    group.by = "CellType",
    split.by = "Condition",
    pt.size  = 0,
    cols     = cond_colors
  ) +
    labs(
      title    = sprintf("%s Expression by Cell Type — All Conditions", gene),
      subtitle = sprintf("All samples | Log-normalised | n = %d cells", ncol(seurat_obj)),
      x        = "Cell Type",
      y        = sprintf("%s expression (log-normalised)", gene)
    ) +
    seurat_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(plot_dir, sprintf("%s_01_Violin_all_conditions.pdf", gene)),
    p_vln_all, width = 18, height = 6
  )

  # --------------------------------------------------------------------------
  # SECTION C: Violin plot — Healing vs Non-Healing, split per cell type
  # --------------------------------------------------------------------------
  cat(sprintf("  [C] VlnPlot Healing vs Non-Healing — %s\n", gene))

  p_vln_hnh <- VlnPlot(
    seurat_hnh,
    features   = gene,
    split.by   = "Condition",
    split.plot = TRUE,
    pt.size    = 0,
    cols       = c("Healing" = "#2196F3", "Non-Healing" = "#F44336")
  ) +
    labs(
      title    = sprintf("%s — Healing vs Non-Healing per Cell Type", gene),
      subtitle = sprintf("Healing: %d cells | Non-Healing: %d cells | Log-normalised",
                         sum(seurat_hnh$Condition == "Healing"),
                         sum(seurat_hnh$Condition == "Non-Healing")),
      x        = "Cell Type",
      y        = sprintf("%s expression (log-normalised)", gene)
    ) +
    seurat_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(plot_dir, sprintf("%s_02_Violin_HvNH.pdf", gene)),
    p_vln_hnh, width = 16, height = 6
  )

  # --------------------------------------------------------------------------
  # SECTION D: Violin plot — solo gruppo Diabetes, per cell type
  # --------------------------------------------------------------------------
  cat(sprintf("  [D] VlnPlot Diabetes group — %s\n", gene))

  p_vln_diab <- VlnPlot(
    seurat_diab,
    features = gene,
    pt.size  = 0.3,
    cols     = ct_colors
  ) +
    labs(
      title    = sprintf("%s Expression per Cell Type — Diabetes Group", gene),
      subtitle = sprintf("n = %d cells | %d diabetic samples | Log-normalised",
                         ncol(seurat_diab),
                         length(unique(seurat_diab$orig.ident))),
      x        = "Cell Type",
      y        = sprintf("%s expression (log-normalised)", gene)
    ) +
    seurat_theme() +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  ggsave(
    file.path(plot_dir, sprintf("%s_03_Violin_Diabetes.pdf", gene)),
    p_vln_diab, width = 14, height = 6
  )

  # --------------------------------------------------------------------------
  # SECTION E: Panel combinato — tutte e tre le figure in colonna
  # --------------------------------------------------------------------------
  cat(sprintf("  [E] Combined panel — %s\n", gene))

  p_panel <- p_vln_all / p_vln_hnh / p_vln_diab +
    plot_annotation(
      title    = sprintf("%s Violin Plot Summary — GSE165816 DFU", gene),
      subtitle = "Top: all conditions | Middle: Healing vs Non-Healing | Bottom: Diabetes group",
      theme    = theme(
        plot.title    = element_text(face = "bold", color = "black", size = 14),
        plot.subtitle = element_text(color = "grey35", size = 11)
      )
    )

  ggsave(
    file.path(plot_dir, sprintf("%s_04_Violin_panel_combined.pdf", gene)),
    p_panel, width = 18, height = 20
  )

  cat(sprintf("  --> Plots salvati per %s\n", gene))
}

# ==============================================================================
# RIEPILOGO FINALE
# ==============================================================================
cat("\n\n=== CCR VIOLIN PLOTS — COMPLETE ===\n")
cat(sprintf("Plots salvati in: %s\n", plot_dir))
cat("Per ogni recettore (CCR1, CCR5, CCR8):\n")
cat("  [GENE]_01_Violin_all_conditions.pdf  — Tutte e tre le condizioni, split per cell type\n")
cat("  [GENE]_02_Violin_HvNH.pdf            — Healing vs Non-Healing, split per cell type\n")
cat("  [GENE]_03_Violin_Diabetes.pdf        — Solo gruppo Diabetes, per cell type\n")
cat("  [GENE]_04_Violin_panel_combined.pdf  — Panel combinato (tutti e tre in colonna)\n")
