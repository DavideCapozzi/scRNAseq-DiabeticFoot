# ==============================================================================
# SCRIPT 03: CCL4 EXPRESSION ANALYSIS
# Dataset: GSE165816 - DFU scRNA-seq
#
# SECTIONS:
#   A. CCL4 expression on global UMAP
#   B. CCL4 expression split by condition (Healing | Non-Healing | Diabetes)
#   C. Differential expression: Healing vs Non-Healing
#      C1. Violin plots per cell type
#      C2. Dot plot: CCL4-CCR5 axis
#      C3. FindMarkers — global
#      C4. FindMarkers — per cell type
#      C5. Pseudobulk DE (statistically correct approach)
#   D. Diabetes control — CCL4 expression in diabetic samples only
#      D1. UMAP feature plot (Diabetes subset)
#      D2. Violin plot per cell type (Diabetes)
#      D3. Three-group comparison box plots
#      D4. Statistical tests (Kruskal-Wallis + pairwise Wilcoxon)
#   E. Summary figure
#
# PREREQUISITE: Script 02 completed
# ==============================================================================

rm(list = ls())
set.seed(42)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(scales)
library(ggrepel)

# ==============================================================================
# 0. SETUP
# ==============================================================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/nature_paper_validation"
res_dir  <- file.path(base_dir, "res")
plot_dir <- file.path(res_dir, "plots", "CCL4")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# GLOBAL PLOT THEME — black axes, no grid, English labels
# ==============================================================================
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
# 1. LOAD ANNOTATED OBJECT
# ==============================================================================
cat("--- Loading annotated object ---\n")
seurat_obj <- readRDS(file.path(res_dir, "seurat_obj_annotated.rds"))

cat(sprintf("Total cells: %d | Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))
cat("Cells per condition:\n"); print(table(seurat_obj$Condition))

n_ct <- length(unique(seurat_obj$CellType))
ct_colors <- scales::hue_pal()(n_ct)
names(ct_colors) <- sort(unique(seurat_obj$CellType))

# ==============================================================================
# CCL4 PRESENCE CHECK
# ==============================================================================
if (!"CCL4" %in% rownames(seurat_obj)) {
  stop("[ERROR] CCL4 not found in features!\n",
       "Check gene names — possible lowercase variant: grep('CCL4', rownames(seurat_obj), value=TRUE)")
}

# Baseline statistics
ccl4_expr <- FetchData(seurat_obj, vars = c("CCL4","Condition","CellType","orig.ident"))
cat("\n=== CCL4 Baseline Statistics ===\n")
cat(sprintf("  CCL4+ cells (>0): %d / %d (%.1f%%)\n",
            sum(ccl4_expr$CCL4 > 0), nrow(ccl4_expr), 100 * mean(ccl4_expr$CCL4 > 0)))
cat("  Mean expression per condition:\n")
print(round(tapply(ccl4_expr$CCL4, ccl4_expr$Condition, mean), 4))
cat("  % CCL4+ cells per condition:\n")
print(round(tapply(ccl4_expr$CCL4 > 0, ccl4_expr$Condition, mean) * 100, 2))

# ==============================================================================
# SECTION A: CCL4 EXPRESSION ON GLOBAL UMAP
# ==============================================================================
cat("\n--- Section A: CCL4 Feature Plot on global UMAP ---\n")
Idents(seurat_obj) <- "CellType"

# A1. CCL4 feature plot (all samples)
p_ccl4_umap <- FeaturePlot(
  seurat_obj, features = "CCL4",
  pt.size = 0.4,
  cols    = c("lightgrey", "#FF9800", "#D62839"),
  order   = TRUE
) +
  labs(title = "CCL4 Expression — Global UMAP",
       subtitle = "All samples | Log-normalised expression",
       x = "UMAP 1", y = "UMAP 2") +
  seurat_theme()

ggsave(file.path(plot_dir, "01_CCL4_FeaturePlot_global.pdf"),
       p_ccl4_umap, width = 8, height = 7)

# A2. Cell type background + CCL4 overlay (side by side)
p_umap_ct_bg <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "CellType",
  label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.3, alpha = 0.4
) +
  scale_color_manual(values = ct_colors) +
  labs(title = "Cell Type Annotation", x = "UMAP 1", y = "UMAP 2") +
  NoLegend() + seurat_theme()

p_ccl4_annotated <- p_umap_ct_bg | p_ccl4_umap
ggsave(file.path(plot_dir, "02_CCL4_UMAP_with_annotation.pdf"),
       p_ccl4_annotated, width = 16, height = 7)

# ==============================================================================
# SECTION B: CCL4 SPLIT BY CONDITION
# ==============================================================================
cat("\n--- Section B: CCL4 split by condition ---\n")

# B1. FeaturePlot split (3 panels, same colour scale)
p_ccl4_split <- FeaturePlot(
  seurat_obj, features = "CCL4",
  split.by   = "Condition",
  pt.size    = 0.3,
  cols       = c("lightgrey", "#D62839"),
  order      = TRUE,
  keep.scale = "all"
) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  plot_annotation(
    title    = "CCL4 Expression by Condition",
    subtitle = "Identical colour scale across panels for direct comparison"
  )
ggsave(file.path(plot_dir, "03_CCL4_FeaturePlot_split_condition.pdf"),
       p_ccl4_split, width = 18, height = 6)

# B2. Violin plot — CCL4 per cell type, split by condition
p_ccl4_vln_all <- VlnPlot(
  seurat_obj, features = "CCL4",
  group.by = "CellType", split.by = "Condition",
  pt.size = 0, cols = cond_colors
) +
  labs(title = "CCL4 Expression by Cell Type — All Conditions",
       x = "Cell Type", y = "CCL4 expression (log-normalised)") +
  seurat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "04_CCL4_Violin_all_conditions.pdf"),
       p_ccl4_vln_all, width = 18, height = 6)

# ==============================================================================
# SECTION C: DIFFERENTIAL EXPRESSION — HEALING vs NON-HEALING
# ==============================================================================
cat("\n--- Section C: CCL4 DE — Healing vs Non-Healing ---\n")

# Subset: exclude Diabetes for this comparison
seurat_hnh <- subset(seurat_obj, subset = Condition %in% c("Healing","Non-Healing"))
cat(sprintf("  Healing: %d cells | Non-Healing: %d cells\n",
            sum(seurat_hnh$Condition == "Healing"),
            sum(seurat_hnh$Condition == "Non-Healing")))

genes_ccl4_axis <- c("CCL4","CCL4L1","CCL4L2","CCR5","CCR1","CCR8")
genes_axis_ok   <- genes_ccl4_axis[genes_ccl4_axis %in% rownames(seurat_obj)]
cat(sprintf("  CCL4 axis genes present: %s\n", paste(genes_axis_ok, collapse = ", ")))

# C1. Violin plot: CCL4 Healing vs Non-Healing per cell type
Idents(seurat_hnh) <- "CellType"
p_ccl4_hnh <- VlnPlot(
  seurat_hnh, features = "CCL4",
  split.by  = "Condition", split.plot = TRUE,
  pt.size   = 0,
  cols      = c("Healing" = "#2196F3", "Non-Healing" = "#F44336")
) +
  labs(title = "CCL4 — Healing vs Non-Healing per Cell Type",
       x = "Cell Type", y = "CCL4 expression (log-normalised)") +
  seurat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "05_CCL4_Violin_HvNH.pdf"), p_ccl4_hnh, width = 16, height = 6)

# C2. Dot plot: CCL4-CCR5 axis
p_dot_hnh <- DotPlot(
  seurat_hnh, features = genes_axis_ok, group.by = "CellType",
  split.by = "Condition",
  cols     = c("Healing" = "#2196F3", "Non-Healing" = "#F44336"),
  dot.scale = 7
) +
  coord_flip() +
  labs(title = "CCL4–CCR5 Axis | Healing vs Non-Healing",
       x = "Gene", y = "Cell Type") +
  seurat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "06_CCL4_axis_DotPlot_HvNH.pdf"), p_dot_hnh, width = 18, height = 6)

# C3. FindMarkers — GLOBAL (Healing vs Non-Healing)
Idents(seurat_hnh) <- "Condition"
deg_ccl4_global <- FindMarkers(
  seurat_hnh,
  ident.1         = "Healing",
  ident.2         = "Non-Healing",
  features        = genes_axis_ok,
  logfc.threshold = 0, min.pct = 0, test.use = "wilcox" #did not put thresholds here
)
deg_ccl4_global$gene <- rownames(deg_ccl4_global)

cat("\n=== CCL4 DE — Global (Healing vs Non-Healing) ===\n")
cat("  avg_log2FC > 0  =  upregulated in Healing\n")
cat("  avg_log2FC < 0  =  downregulated in Healing\n\n")
print(deg_ccl4_global[, c("gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj")])

if ("CCL4" %in% rownames(deg_ccl4_global)) {
  r   <- deg_ccl4_global["CCL4", ]
  dir <- if (r$avg_log2FC > 0) "UPREGULATED in Healing" else "DOWNREGULATED in Healing"
  sig <- if (r$p_val_adj < 0.05) "(FDR < 0.05 — SIGNIFICANT)" else "(not significant)"
  cat(sprintf("\n>>> CCL4: %s | log2FC = %.3f | FDR = %.4f %s\n",
              dir, r$avg_log2FC, r$p_val_adj, sig))
}
write.csv(deg_ccl4_global, file.path(res_dir, "CCL4_DE_global_HvNH.csv"), row.names = FALSE)

# C4. FindMarkers — per cell type
cat("\n--- CCL4 DE per cell type ---\n")
ccl4_de_per_ct <- list()

for (ct in unique(seurat_hnh$CellType)) {
  sub  <- subset(seurat_hnh, subset = CellType == ct)
  n_h  <- sum(sub$Condition == "Healing")
  n_nh <- sum(sub$Condition == "Non-Healing")
  if (n_h < 10 || n_nh < 10) {
    cat(sprintf("  [SKIP] %s: insufficient cells (H=%d, NH=%d)\n", ct, n_h, n_nh)); next
  }
  Idents(sub) <- "Condition"
  tryCatch({
    deg <- FindMarkers(sub, ident.1 = "Healing", ident.2 = "Non-Healing",
                       features = genes_axis_ok, logfc.threshold = 0,
                       min.pct = 0, test.use = "wilcox")
    deg$gene <- rownames(deg); deg$CellType <- ct
    ccl4_de_per_ct[[ct]] <- deg
    if ("CCL4" %in% rownames(deg)) {
      r   <- deg["CCL4", ]
      dir <- if (r$avg_log2FC > 0) "UP  " else "DOWN"
      sig <- if (r$p_val_adj < 0.05) "**" else "ns"
      cat(sprintf("  %s | %s log2FC=%.3f %s\n", ct, dir, r$avg_log2FC, sig))
    }
  }, error = function(e) cat(sprintf("  [ERROR] %s: %s\n", ct, e$message)))
}

if (length(ccl4_de_per_ct) > 0) {
  deg_all_ct <- bind_rows(ccl4_de_per_ct)
  write.csv(deg_all_ct, file.path(res_dir, "CCL4_DE_per_CellType.csv"), row.names = FALSE)
  cat("\n=== CCL4 significant per cell type (FDR < 0.05) ===\n")
  sig_ct <- deg_all_ct %>%
    filter(gene == "CCL4", p_val_adj < 0.05) %>%
    mutate(Direction = ifelse(avg_log2FC > 0, "UP in Healing", "DOWN in Healing")) %>%
    select(CellType, gene, Direction, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    arrange(p_val_adj)
  print(sig_ct)
}

# C5. PSEUDOBULK DE — statistically correct approach
# Aggregates cells per patient per cell type → avoids pseudoreplication
cat("\n--- Pseudobulk DE — CCL4 (Healing vs Non-Healing) ---\n")
cat("  Aggregation by patient removes pseudoreplication bias\n")

pseudobulk_results <- list()

for (ct in unique(seurat_hnh$CellType)) {
  ct_cells <- subset(seurat_hnh, subset = CellType == ct)
  pb_mat   <- tryCatch(
    AggregateExpression(ct_cells, group.by = "orig.ident",
                        slot = "counts", return.seurat = FALSE)$RNA,
    error = function(e) NULL
  )
  if (is.null(pb_mat) || ncol(pb_mat) < 4 || !"CCL4" %in% rownames(pb_mat)) next

  patient_cond <- seurat_hnh@meta.data %>% select(orig.ident, Condition) %>% distinct()
  cols_ok      <- intersect(colnames(pb_mat), patient_cond$orig.ident)
  if (length(cols_ok) < 4) next

  pb_sub   <- pb_mat[, cols_ok, drop = FALSE]
  cond_sub <- patient_cond$Condition[match(cols_ok, patient_cond$orig.ident)]

  # log-CPM normalisation
  cpm_ccl4 <- log1p(pb_sub["CCL4", ] / (colSums(pb_sub) / 1e6))
  h_vals   <- cpm_ccl4[cond_sub == "Healing"]
  nh_vals  <- cpm_ccl4[cond_sub == "Non-Healing"]

  if (length(h_vals) >= 2 && length(nh_vals) >= 2) {
    tt <- t.test(h_vals, nh_vals)
    pseudobulk_results[[ct]] <- data.frame(
      CellType     = ct,
      mean_Healing = round(mean(h_vals), 4),
      mean_NonHeal = round(mean(nh_vals), 4),
      log2FC       = round(mean(h_vals) - mean(nh_vals), 4),
      pvalue       = tt$p.value,
      n_Healing    = length(h_vals),
      n_NonHealing = length(nh_vals)
    )
  }
}

if (length(pseudobulk_results) > 0) {
  pb_df      <- bind_rows(pseudobulk_results)
  pb_df$padj <- p.adjust(pb_df$pvalue, method = "BH")
  pb_df      <- pb_df[order(pb_df$pvalue), ]

  cat("\n=== PSEUDOBULK CCL4 — Results ===\n")
  cat("  log2FC > 0 = upregulated in Healing\n\n")
  print(pb_df)
  write.csv(pb_df, file.path(res_dir, "CCL4_pseudobulk_HvNH.csv"), row.names = FALSE)

  p_pb <- ggplot(pb_df, aes(x = reorder(CellType, log2FC), y = log2FC,
                              fill = ifelse(padj < 0.05, "FDR < 0.05", "n.s."))) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("FDR < 0.05" = "#D62839", "n.s." = "grey65"),
                       name = "Significance") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    coord_flip() +
    labs(title = "CCL4 Pseudobulk DE: Healing vs Non-Healing",
         subtitle = "log2FC > 0 = upregulated in Healing | red = FDR < 0.05",
         x = "Cell Type", y = "log2 Fold Change (pseudobulk log-CPM)") +
    theme_dfu(12)
  ggsave(file.path(plot_dir, "07_CCL4_pseudobulk_barplot.pdf"), p_pb, width = 10, height = 7)
}

# ==============================================================================
# SECTION D: DIABETES CONTROL
# ==============================================================================
cat("\n--- Section D: CCL4 in Diabetes group (control) ---\n")

seurat_diab <- subset(seurat_obj, subset = Condition == "Diabetes")
cat(sprintf("  Diabetes cells: %d | Samples: %s\n",
            ncol(seurat_diab),
            paste(unique(seurat_diab$orig.ident), collapse = ", ")))

# D1. FeaturePlot — Diabetes subset UMAP
p_ccl4_diab_umap <- FeaturePlot(
  seurat_diab, features = "CCL4",
  pt.size = 0.5,
  cols    = c("lightgrey", "#FF9800", "#D62839"),
  order   = TRUE
) +
  labs(title = "CCL4 Expression — Diabetes Group (UMAP)",
       subtitle = sprintf("Diabetic control | No active ulcer | n = %d cells",
                          ncol(seurat_diab)),
       x = "UMAP 1", y = "UMAP 2") +
  seurat_theme()
ggsave(file.path(plot_dir, "08_CCL4_Diabetes_FeaturePlot.pdf"),
       p_ccl4_diab_umap, width = 8, height = 7)

# D2. Violin plot — CCL4 per cell type in Diabetes
Idents(seurat_diab) <- "CellType"
p_ccl4_diab_vln <- VlnPlot(
  seurat_diab, features = "CCL4",
  pt.size = 0.3, cols = ct_colors
) +
  labs(title = "CCL4 Expression per Cell Type — Diabetes Group",
       subtitle = sprintf("n = %d cells | %d diabetic samples",
                          ncol(seurat_diab), length(unique(seurat_diab$orig.ident))),
       x = "Cell Type", y = "CCL4 expression (log-normalised)") +
  seurat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(plot_dir, "09_CCL4_Diabetes_Violin.pdf"),
       p_ccl4_diab_vln, width = 14, height = 6)

# D3. Three-group comparison
cat("\n--- Three-group comparison: Healing / Non-Healing / Diabetes ---\n")

ccl4_3grp <- FetchData(seurat_obj,
                        vars = c("CCL4","Condition","CellType","orig.ident"))

# Aggregate per patient (removes pseudoreplication for boxplots)
ccl4_patient <- ccl4_3grp %>%
  group_by(orig.ident, Condition) %>%
  summarise(mean_CCL4 = mean(CCL4),
            pct_pos   = mean(CCL4 > 0) * 100,
            .groups = "drop")

# Box plot: mean CCL4 per patient per condition
p_box_global <- ggplot(ccl4_patient, aes(x = Condition, y = mean_CCL4, fill = Condition)) +
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.1, size = 2.5, color = "black", alpha = 0.7) +
  scale_fill_manual(values = cond_colors) +
  labs(title = "CCL4 Mean Expression — Three Conditions",
       subtitle = "One dot = one patient | log-normalised",
       x = "Condition", y = "Mean CCL4 expression") +
  theme_dfu(13) +
  theme(legend.position = "none")

# Box plot: % CCL4+ cells per patient
p_box_pct <- ggplot(ccl4_patient, aes(x = Condition, y = pct_pos, fill = Condition)) +
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.1, size = 2.5, color = "black", alpha = 0.7) +
  scale_fill_manual(values = cond_colors) +
  labs(title = "% CCL4+ Cells — Three Conditions",
       subtitle = "One dot = one patient",
       x = "Condition", y = "% CCL4-positive cells") +
  theme_dfu(13) +
  theme(legend.position = "none")

p_box_panel <- p_box_global | p_box_pct
ggsave(file.path(plot_dir, "10_CCL4_ThreeGroups_boxplot.pdf"),
       p_box_panel, width = 12, height = 6)

# Box plot: CCL4 per cell type — all three conditions
ccl4_ct_patient <- ccl4_3grp %>%
  group_by(orig.ident, Condition, CellType) %>%
  summarise(mean_CCL4 = mean(CCL4), .groups = "drop")

p_box_celltype <- ggplot(ccl4_ct_patient,
                          aes(x = CellType, y = mean_CCL4, fill = Condition)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(width = 0.75),
               color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cond_colors) +
  labs(title = "CCL4 Mean Expression per Cell Type — Three Conditions",
       subtitle = "Aggregated per patient",
       x = "Cell Type", y = "Mean CCL4 (log-normalised)") +
  theme_dfu(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "11_CCL4_CellType_ThreeGroups_boxplot.pdf"),
       p_box_celltype, width = 18, height = 7)

# D4. Statistical tests — three-group comparison
cat("\n=== Statistical tests: CCL4 across three conditions ===\n")

kw_test <- kruskal.test(mean_CCL4 ~ Condition, data = ccl4_patient)
cat(sprintf("Kruskal-Wallis: chi2 = %.3f, df = %d, p = %.4f\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

pw_test <- pairwise.wilcox.test(ccl4_patient$mean_CCL4,
                                 ccl4_patient$Condition,
                                 p.adjust.method = "BH")
cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
print(pw_test$p.value)

# ==============================================================================
# SECTION E: SUMMARY FIGURE
# ==============================================================================
cat("\n--- Section E: Summary figure ---\n")

p_summary <- (p_ccl4_umap | p_ccl4_split) /
             (p_box_global | p_box_pct) +
  plot_annotation(
    title    = "CCL4 Analysis Summary — GSE165816 DFU",
    subtitle = "Top: global UMAP expression (left) and split by condition (right) | Bottom: patient-level statistics"
  )
ggsave(file.path(plot_dir, "12_CCL4_summary_figure.pdf"),
       p_summary, width = 20, height = 14)

# ==============================================================================
# EXPORT RESULTS
# ==============================================================================
cat("\n--- Exporting summary tables ---\n")

# Summary: mean, %, per cell type per condition
summary_table <- ccl4_3grp %>%
  group_by(Condition, CellType) %>%
  summarise(
    n_cells      = n(),
    mean_expr    = round(mean(CCL4), 4),
    median_expr  = round(median(CCL4), 4),
    pct_positive = round(mean(CCL4 > 0) * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(CellType, Condition)

write.csv(summary_table,
          file.path(res_dir, "CCL4_summary_CellType_Condition.csv"),
          row.names = FALSE)

write.csv(ccl4_patient,
          file.path(res_dir, "CCL4_per_patient_3groups.csv"),
          row.names = FALSE)

cat("\n=== CCL4 ANALYSIS COMPLETE ===\n")
cat("Results saved:\n")
cat("  res/CCL4_DE_global_HvNH.csv            — Global DE Healing vs Non-Healing\n")
cat("  res/CCL4_DE_per_CellType.csv           — DE per cell type\n")
cat("  res/CCL4_pseudobulk_HvNH.csv           — Pseudobulk DE\n")
cat("  res/CCL4_summary_CellType_Condition.csv — Statistics per cell type\n")
cat("  res/CCL4_per_patient_3groups.csv        — Patient-level means (3 groups)\n")
cat("\nPlots saved in:", plot_dir, "\n")
cat("  01 — CCL4 global UMAP\n")
cat("  02 — CCL4 UMAP with annotation\n")
cat("  03 — CCL4 split by condition\n")
cat("  04 — Violin: all conditions\n")
cat("  05 — Violin: Healing vs Non-Healing\n")
cat("  06 — Dot plot: CCL4 axis\n")
cat("  07 — Pseudobulk bar plot\n")
cat("  08 — Diabetes FeaturePlot\n")
cat("  09 — Diabetes Violin\n")
cat("  10 — Three-group box plots (patient-level)\n")
cat("  11 — Three-group box plots (per cell type)\n")
cat("  12 — Summary figure\n")
