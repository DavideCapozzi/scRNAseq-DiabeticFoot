# ==============================================================================
# CCR5 Expression Level Assessment
# Purpose: Determine whether CCR5 can be described as "highly expressed"
#          by comparing its expression to the global gene expression distribution.
# Output:  Summary stats table, per-cell-type table, plots, and a plain-text
#          verdict suitable for paper writing.
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

# ------------------------------------------------------------------------------
# PATHS  (Windows / Davide's PC)
# ------------------------------------------------------------------------------
seurat_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds"

output_dir  <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/res_validation/CCR5_expression_assessment/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

target_gene <- "CCR5"

# ------------------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------------------
message("Loading Seurat object ...")
#seurat_obj <- readRDS(seurat_path)
DefaultAssay(seurat_obj) <- "RNA"

if (!target_gene %in% rownames(seurat_obj)) {
  stop(paste(target_gene, "not found in the RNA assay. Check gene name capitalisation."))
}

# Make sure cell-type labels are set as identity
# (the annotated object should have a column called "celltype" or similar)
ct_col <- intersect(c("celltype", "cell_type", "CellType", "manual_annotation",
                       "celltype_id", "celltypeID"),
                    colnames(seurat_obj@meta.data))[1]

if (!is.na(ct_col)) {
  Idents(seurat_obj) <- ct_col
  message("Cell-type identity set to column: ", ct_col)
} else {
  message("No recognised cell-type column found – using existing Idents()")
}

# ------------------------------------------------------------------------------
# 2. EXTRACT NORMALISED EXPRESSION MATRIX
# ------------------------------------------------------------------------------
message("Extracting normalised expression data ...")

# Use the 'data' slot (log-normalised counts)
expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# CCR5 vector across all cells
CCR5_expr <- as.numeric(expr_matrix[target_gene, ])

# ------------------------------------------------------------------------------
# 3. GLOBAL GENE-LEVEL STATS (compare CCR5 against ALL genes)
# ------------------------------------------------------------------------------
message("Computing gene-level summary statistics across all ", ncol(expr_matrix), " cells ...")

# For each gene: mean expression and fraction of cells with detectable expression
# We do this in chunks to avoid RAM issues on 292 k cells
chunk_size  <- 5000
gene_names  <- rownames(expr_matrix)
n_genes     <- length(gene_names)
n_chunks    <- ceiling(n_genes / chunk_size)

gene_means  <- numeric(n_genes)
gene_pct    <- numeric(n_genes)

for (i in seq_len(n_chunks)) {
  idx_start <- (i - 1) * chunk_size + 1
  idx_end   <- min(i * chunk_size, n_genes)
  chunk     <- expr_matrix[idx_start:idx_end, , drop = FALSE]
  gene_means[idx_start:idx_end] <- Matrix::rowMeans(chunk)
  gene_pct  [idx_start:idx_end] <- Matrix::rowMeans(chunk > 0)
}
names(gene_means) <- gene_names
names(gene_pct)   <- gene_names

# CCR5-specific values
CCR5_mean <- gene_means[target_gene]
CCR5_pct  <- gene_pct [target_gene]

# Ranks (higher rank = more expressed)
rank_by_mean <- rank(gene_means, ties.method = "average")
rank_by_pct  <- rank(gene_pct,   ties.method = "average")

CCR5_mean_rank       <- rank_by_mean[target_gene]
CCR5_mean_percentile <- CCR5_mean_rank / n_genes * 100

CCR5_pct_rank        <- rank_by_pct[target_gene]
CCR5_pct_percentile  <- CCR5_pct_rank / n_genes * 100

# Global summary
global_summary <- data.frame(
  Metric        = c("Mean normalised expression",
                    "% cells expressing CCR5 (> 0)",
                    "Rank by mean expression (out of all genes)",
                    "Percentile by mean expression",
                    "Rank by % expressing cells",
                    "Percentile by % expressing cells",
                    "Total genes in dataset",
                    "Total cells in dataset"),
  CCR5_value    = c(round(CCR5_mean, 4),
                    round(CCR5_pct * 100, 2),
                    round(CCR5_mean_rank),
                    round(CCR5_mean_percentile, 1),
                    round(CCR5_pct_rank),
                    round(CCR5_pct_percentile, 1),
                    n_genes,
                    ncol(expr_matrix))
)

print(global_summary)

# ------------------------------------------------------------------------------
# 4. PER-CELL-TYPE STATS FOR CCR5
# ------------------------------------------------------------------------------
message("Computing per-cell-type stats for CCR5 ...")

cell_types <- as.character(Idents(seurat_obj))

ct_stats <- data.frame(
  cell_type  = cell_types,
  expression = CCR5_expr
) %>%
  group_by(cell_type) %>%
  summarise(
    n_cells          = n(),
    mean_expr        = round(mean(expression), 4),
    median_expr      = round(median(expression), 4),
    pct_expressing   = round(mean(expression > 0) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_expr))

print(ct_stats)

# ------------------------------------------------------------------------------
# 5. REFERENCE COMPARISON: HIGH vs LOW expressed genes
# ------------------------------------------------------------------------------
# Identify top-50 expressed genes globally as "housekeeping / highly expressed"
# and bottom-50 as poorly expressed
message("Identifying reference gene sets ...")

sorted_genes    <- sort(gene_means, decreasing = TRUE)
top50_genes     <- names(sorted_genes)[1:50]
bottom50_genes  <- names(sorted_genes)[(n_genes - 49):n_genes]

top50_mean_avg    <- mean(gene_means[top50_genes])
bottom50_mean_avg <- mean(gene_means[bottom50_genes])

# Z-score of CCR5 against global distribution
z_score_CCR5 <- (CCR5_mean - mean(gene_means)) / sd(gene_means)

reference_summary <- data.frame(
  Reference_set            = c("CCR5",
                                "Average of top-50 expressed genes",
                                "Average of bottom-50 expressed genes",
                                "Global gene mean expression"),
  Mean_normalised_expr     = c(round(CCR5_mean, 4),
                                round(top50_mean_avg, 4),
                                round(bottom50_mean_avg, 6),
                                round(mean(gene_means), 4))
)

print(reference_summary)

# ------------------------------------------------------------------------------
# 6. PLOTS
# ------------------------------------------------------------------------------
message("Generating plots ...")

## 6a. Density plot of gene means – CCR5 highlighted
gene_df <- data.frame(gene = gene_names, mean_expr = gene_means)

p_density <- ggplot(gene_df, aes(x = mean_expr)) +
  geom_density(fill = "steelblue", alpha = 0.4, colour = "steelblue") +
  geom_vline(xintercept = CCR5_mean, colour = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = CCR5_mean, y = Inf,
           label = paste0(" CCR5\n(mean = ", round(CCR5_mean, 3), ")"),
           colour = "red", hjust = 0, vjust = 1.5, size = 3.5) +
  scale_x_log10() +
  labs(title = "Distribution of mean normalised expression across all genes",
       subtitle = paste0("CCR5 is at the ", round(CCR5_mean_percentile, 1),
                         "th percentile (n = ", n_genes, " genes, ",
                         ncol(expr_matrix), " cells)"),
       x = "Mean normalised expression (log10 scale)",
       y = "Density") +
  theme_classic(base_size = 13)

ggsave(file.path(output_dir, "1_CCR5_vs_global_gene_expression_density.png"),
       p_density, width = 9, height = 5, dpi = 300)

## 6b. Barplot of CCR5 mean expression per cell type
p_bar <- ggplot(ct_stats, aes(x = reorder(cell_type, mean_expr), y = mean_expr,
                               fill = pct_expressing)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(pct_expressing, "%")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_gradient(low = "lightyellow", high = "darkred",
                      name = "% cells\nexpressing") +
  labs(title = "CCR5 mean normalised expression by cell type",
       subtitle = "Labels = % cells with detectable expression",
       x = NULL, y = "Mean normalised expression") +
  theme_classic(base_size = 12) +
  expand_limits(y = max(ct_stats$mean_expr) * 1.15)

ggsave(file.path(output_dir, "2_CCR5_mean_expression_per_celltype.png"),
       p_bar, width = 10, height = 7, dpi = 300)

## 6c. Violin plot of CCR5 expression per cell type
p_vln <- VlnPlot(seurat_obj, features = target_gene, pt.size = 0,
                  assay = "RNA", slot = "data") +
  labs(title = paste0("CCR5 expression per cell type"),
       y = "Normalised expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "3_CCR5_violin_per_celltype.png"),
       p_vln, width = 12, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 7. SAVE TABLES TO EXCEL
# ------------------------------------------------------------------------------
message("Saving Excel summary ...")

wb <- createWorkbook()

addWorksheet(wb, "Global_summary")
writeData(wb, "Global_summary", global_summary)

addWorksheet(wb, "Per_celltype_stats")
writeData(wb, "Per_celltype_stats", ct_stats)

addWorksheet(wb, "Reference_comparison")
writeData(wb, "Reference_comparison", reference_summary)

addWorksheet(wb, "Top50_genes")
writeData(wb, "Top50_genes", data.frame(Gene = top50_genes,
                                         Mean_expr = round(gene_means[top50_genes], 4)))

saveWorkbook(wb, file.path(output_dir, "CCR5_expression_assessment.xlsx"),
             overwrite = TRUE)

# ------------------------------------------------------------------------------
# 8. PLAIN-TEXT VERDICT
# ------------------------------------------------------------------------------
verdict_file <- file.path(output_dir, "CCR5_expression_verdict.txt")

verdict_text <- paste0(
  "==========================================================\n",
  "  CCR5 EXPRESSION LEVEL ASSESSMENT – VERDICT\n",
  "  Generated: ", Sys.time(), "\n",
  "==========================================================\n\n",
  "Dataset:  ", ncol(expr_matrix), " cells  |  ", n_genes, " genes\n\n",
  "--- CCR5 global stats ---\n",
  "  Mean normalised expression : ", round(CCR5_mean, 4), "\n",
  "  % cells expressing (> 0)  : ", round(CCR5_pct * 100, 2), "%\n",
  "  Percentile (by mean expr)  : ", round(CCR5_mean_percentile, 1), "th\n",
  "  Z-score vs all genes       : ", round(z_score_CCR5, 2), "\n\n",
  "--- Reference ---\n",
  "  Avg mean expr of TOP-50 genes     : ", round(top50_mean_avg, 4), "\n",
  "  Avg mean expr of BOTTOM-50 genes  : ", round(bottom50_mean_avg, 6), "\n",
  "  Global mean expr (all genes)      : ", round(mean(gene_means), 4), "\n\n",
  "--- Cell type with highest CCR5 expression ---\n",
  "  ", ct_stats$cell_type[1], "\n",
  "  Mean expr = ", ct_stats$mean_expr[1],
  " | % expressing = ", ct_stats$pct_expressing[1], "%\n\n",
  "--- VERDICT ---\n",
  if (CCR5_mean_percentile >= 75) {
    paste0("CCR5 is in the TOP ", round(100 - CCR5_mean_percentile, 1),
           "% of expressed genes (", round(CCR5_mean_percentile, 1),
           "th percentile).\n",
           "--> CCR5 can be described as HIGHLY EXPRESSED in the dataset.\n")
  } else if (CCR5_mean_percentile >= 50) {
    paste0("CCR5 is above median expression (", round(CCR5_mean_percentile, 1),
           "th percentile) but not in the top 25%.\n",
           "--> Consider describing CCR5 as MODERATELY EXPRESSED, or\n",
           "    highlight high expression specifically in the cell type\n",
           "    where it peaks (see per-cell-type table).\n")
  } else {
    paste0("CCR5 is below median expression globally (",
           round(CCR5_mean_percentile, 1), "th percentile).\n",
           "--> Globally CCR5 is LOWLY EXPRESSED, but check whether it is\n",
           "    highly expressed within a specific cell type before claiming\n",
           "    high expression in the paper.\n")
  },
  "\n==========================================================\n"
)

writeLines(verdict_text, verdict_file)
message(verdict_text)

message("\nAll outputs saved to: ", output_dir)
message("Done.")
