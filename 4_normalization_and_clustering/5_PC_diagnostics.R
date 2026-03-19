# ============================================================
#  PC_diagnostics.R
#  Diagnostics per la selezione del numero ottimale di PC
#  Oggetto: seurat_obj (già caricato in sessione)
#  Output:  ./PC_diagnostics/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)

seurat_file <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/5_annotation/manual_annotation_res/seurat_res_0.7_with_celltypeID.rds"
seurat_obj <- readRDS(seurat_file)

# ── output dir ───────────────────────────────────────────────
out_dir <- "G:/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/4_normalization_and_clustering/PC_diagnostics/"
dir.create(out_dir, showWarnings = FALSE)

# ── costanti ─────────────────────────────────────────────────
N_PCS_USED  <- 6          # quello attuale
N_PCS_TOTAL <- 50         # npcs usati in RunPCA
CUTOFFS     <- c(6, 10, 15, 20)   # cutoff da confrontare nei plot


# ============================================================
#  UTILITY: salva ggplot con dimensioni consistenti
# ============================================================
save_plot <- function(p, filename, w = 10, h = 7) {
  ggsave(file.path(out_dir, filename), plot = p,
         width = w, height = h, dpi = 200, bg = "white")
  message("  → salvato: ", filename)
}


# ============================================================
#  1. ELBOW PLOT — varianza % per PC + varianza cumulativa
#     con marcatori per ogni cutoff candidato
# ============================================================
message("[1/6] Elbow plot con varianza cumulativa...")

pca       <- seurat_obj@reductions$pca
stdev     <- pca@stdev
total_var <- pca@misc$total.variance
var_exp   <- (stdev^2) / total_var * 100
cum_var   <- cumsum(var_exp)

df_var <- data.frame(
  PC      = seq_len(N_PCS_TOTAL),
  var_exp = var_exp,
  cum_var = cum_var,
  used    = seq_len(N_PCS_TOTAL) <= N_PCS_USED
)

# calcola il "delta" tra PC consecutive (utile per vedere il gomito)
df_var$delta <- c(NA, -diff(var_exp))

p_elbow <- ggplot(df_var, aes(x = PC, y = var_exp, fill = used)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_text(aes(label = ifelse(PC <= 20, sprintf("%.1f", var_exp), "")),
            vjust = -0.4, size = 2.4, color = "grey25") +
  # linea varianza cumulativa (asse secondario)
  geom_line(aes(x = PC, y = cum_var / 4), color = "grey40",
            linewidth = 0.6, linetype = "solid", inherit.aes = FALSE,
            data = df_var) +
  geom_point(aes(x = PC, y = cum_var / 4), color = "grey40",
             size = 1.2, inherit.aes = FALSE, data = df_var) +
  # marcatori per i cutoff candidati
  geom_vline(xintercept = CUTOFFS + 0.5,
             linetype   = c("dashed", "dotted", "dotdash", "longdash"),
             color      = c("firebrick","#e07b00","#1a7abd","#4caf50"),
             linewidth  = 0.8) +
  annotate("text",
           x     = CUTOFFS + 0.6,
           y     = max(var_exp) * 0.85,
           label = paste0("PC", CUTOFFS,
                          "\n", sprintf("%.1f%%", cum_var[CUTOFFS])),
           hjust = 0, size = 3,
           color = c("firebrick","#e07b00","#1a7abd","#4caf50")) +
  scale_fill_manual(values = c("TRUE" = "#2171b5", "FALSE" = "#bdd7e7")) +
  scale_x_continuous(breaks = c(1, CUTOFFS, seq(10, N_PCS_TOTAL, 10)) |> unique() |> sort()) +
  scale_y_continuous(
    name     = "% Variance explained",
    expand   = expansion(mult = c(0, 0.15)),
    sec.axis = sec_axis(~ . * 4, name = "Cumulative variance (%)")
  ) +
  labs(
    title    = "PCA — Variance explained per PC",
    subtitle = sprintf("PC1–6 = %.1f%%  |  PC1–10 = %.1f%%  |  PC1–20 = %.1f%%  |  All 50 PCs = %.1f%%",
                       cum_var[6], cum_var[10], cum_var[20], cum_var[50]),
    x = "Principal Component"
  ) +
  theme_classic(base_size = 13)

save_plot(p_elbow, "01_elbow_plot.png", w = 13, h = 6)


# ============================================================
#  2. DELTA VARIANCE — derivata prima della varianza spiegata
#     Il gomito vero è dove il delta diventa piatto (< soglia)
# ============================================================
message("[2/6] Delta variance (derivata prima)...")

df_delta <- df_var |> filter(!is.na(delta))

p_delta <- ggplot(df_delta, aes(x = PC, y = delta)) +
  geom_col(aes(fill = delta > 0.2), width = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0.2, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  annotate("text", x = N_PCS_TOTAL - 2, y = 0.25,
           label = "δ = 0.2 threshold", color = "firebrick", size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#aec7e8")) +
  scale_x_continuous(breaks = c(2, CUTOFFS, seq(10, N_PCS_TOTAL, 10)) |> unique() |> sort()) +
  labs(
    title    = "Drop in variance explained between consecutive PCs (Δ)",
    subtitle = "Bars above threshold = PC still captures meaningful signal",
    x = "Principal Component",
    y = "Δ variance explained (%)"
  ) +
  theme_classic(base_size = 13)

save_plot(p_delta, "02_delta_variance.png", w = 13, h = 5)


# ============================================================
#  3. HEATMAP TOP FEATURE LOADINGS — PC1:15
#     Mostra se le PC tardive (7-15) hanno ancora struttura biologica
#     o sono rumore
# ============================================================
message("[3/6] Heatmap feature loadings PC1-15...")

# DimHeatmap è la funzione nativa Seurat per questo scopo
# Salviamo come PNG aprendo device manualmente (DimHeatmap non è ggplot)
png(file.path(out_dir, "03_dimheatmap_PC1-15.png"),
    width = 2400, height = 3200, res = 200)
DimHeatmap(seurat_obj,
           dims       = 1:15,
           cells      = 500,    # subsample per velocità
           balanced   = TRUE,
           fast       = FALSE)
dev.off()
message("  → salvato: 03_dimheatmap_PC1-15.png")


# ============================================================
#  4. PC LOADINGS DOTPLOT — top 10 geni per PC1:10
#     Permette di vedere se le PC 7-10 hanno geni biologicamente
#     sensati o sono noise
# ============================================================
message("[4/6] Top gene loadings per PC1-10...")

# Estrai feature loadings
loadings <- pca@feature.loadings  # matrice [2000 genes x 50 PCs]

top_genes <- lapply(1:10, function(pc) {
  vals <- loadings[, pc]
  top  <- sort(abs(vals), decreasing = TRUE)[1:10]
  data.frame(
    PC        = paste0("PC", pc),
    PC_num    = pc,
    gene      = names(top),
    loading   = vals[names(top)],
    abs_load  = top
  )
}) |> bind_rows()

top_genes$gene <- factor(top_genes$gene,
                         levels = top_genes |>
                           arrange(PC_num, desc(abs_load)) |>
                           pull(gene) |> unique())

p_loadings <- ggplot(top_genes, aes(x = PC, y = gene, size = abs_load,
                                    color = loading)) +
  geom_point() +
  scale_color_gradient2(low = "#2171b5", mid = "grey90", high = "#cb181d",
                        midpoint = 0, name = "Loading") +
  scale_size_continuous(range = c(1, 6), name = "|Loading|") +
  labs(
    title    = "Top 10 gene loadings per PC (PC1–10)",
    subtitle = "PCs with biologically meaningful genes → real signal, not noise",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_loadings, "04_gene_loadings_PC1-10.png", w = 12, h = 14)


# ============================================================
#  5. UMAP SENSITIVITY — confronto UMAP con dims 6 / 10 / 15
#     colorato per seurat_clusters attuali
#     NOTA: ricalcola FindNeighbors + RunUMAP in locale,
#           NON sovrascrive l'oggetto originale
# ============================================================
message("[5/6] UMAP sensitivity analysis (6 / 10 / 15 PCs)...")

umap_plots <- list()

for (nd in c(6, 10, 15)) {
  tmp <- FindNeighbors(seurat_obj, dims = 1:nd, verbose = FALSE)
  tmp <- RunUMAP(tmp, dims = 1:nd, verbose = FALSE,
                 reduction.name = "umap_tmp",
                 reduction.key  = "umaptmp_")
  
  df_umap <- as.data.frame(tmp@reductions$umap_tmp@cell.embeddings)
  df_umap$cluster <- seurat_obj$seurat_clusters
  
  umap_plots[[as.character(nd)]] <- ggplot(df_umap,
                                           aes(x = umaptmp_1, y = umaptmp_2, color = cluster)) +
    geom_point(size = 0.05, alpha = 0.3) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title   = sprintf("UMAP  —  dims = 1:%d", nd),
         subtitle = sprintf("Cumulative variance = %.1f%%", cum_var[nd]),
         x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(legend.text  = element_text(size = 7),
          legend.title = element_blank())
  
  rm(tmp); gc()
}

p_umap_sensitivity <- wrap_plots(umap_plots, ncol = 3) +
  plot_annotation(
    title    = "UMAP Sensitivity Analysis",
    subtitle = "Same clusters (from dims=1:6) proiettati su UMAP ricalcolate con diversi numeri di PC",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11))
  )

save_plot(p_umap_sensitivity, "05_umap_sensitivity.png", w = 18, h = 6)


# ============================================================
#  6. CLUSTER STABILITY TABLE — quante cellule cambiano cluster
#     tra risoluzione calcolata con 6 vs 10 PC?
# ============================================================
message("[6/6] Cluster stability: concordanza 6 vs 10 PC...")

tmp10 <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
tmp10 <- FindClusters(tmp10, resolution = 0.7, verbose = FALSE,
                      cluster.name = "clusters_10pc")

df_comp <- data.frame(
  clusters_6pc  = as.character(seurat_obj$seurat_clusters),
  clusters_10pc = as.character(tmp10$clusters_10pc)
)

# % concordanza globale
concordance <- mean(df_comp$clusters_6pc == df_comp$clusters_10pc)

# contingency heatmap
cont_mat <- df_comp |>
  count(clusters_6pc, clusters_10pc) |>
  group_by(clusters_6pc) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup()

p_stability <- ggplot(cont_mat,
                      aes(x = clusters_10pc, y = clusters_6pc, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct > 2, sprintf("%.0f%%", pct), "")),
            size = 2.8, color = "white") +
  scale_fill_gradient(low = "white", high = "#08306b",
                      name = "% cells\n(row-norm.)") +
  labs(
    title    = "Cluster concordance: 6 PC vs 10 PC",
    subtitle = sprintf("Global concordance = %.1f%%  |  Resolution = 0.7 in entrambi i casi",
                       concordance * 100),
    x = "Clusters (dims = 1:10)",
    y = "Clusters (dims = 1:6)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(size = 9))

save_plot(p_stability, "06_cluster_stability_6vs10PC.png", w = 11, h = 9)

rm(tmp10); gc()


# ============================================================
#  SUMMARY REPORT — stampa a console
# ============================================================
cat("\n")
cat("══════════════════════════════════════════════\n")
cat("  PC DIAGNOSTICS — SUMMARY\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Cells           : %d\n", ncol(seurat_obj)))
cat(sprintf("  Genes           : %d\n", nrow(seurat_obj)))
cat(sprintf("  Variable features: %d\n", length(VariableFeatures(seurat_obj))))
cat(sprintf("  PCs computed    : %d\n", N_PCS_TOTAL))
cat("----------------------------------------------\n")
cat(sprintf("  Variance @ PC6  : %.2f%%\n", cum_var[6]))
cat(sprintf("  Variance @ PC10 : %.2f%%\n", cum_var[10]))
cat(sprintf("  Variance @ PC15 : %.2f%%\n", cum_var[15]))
cat(sprintf("  Variance @ PC20 : %.2f%%\n", cum_var[20]))
cat(sprintf("  Variance @ PC50 : %.2f%%\n", cum_var[50]))
cat("----------------------------------------------\n")
cat(sprintf("  Cluster concordance 6 vs 10 PC: %.1f%%\n", concordance * 100))
cat("══════════════════════════════════════════════\n")
cat("  Output salvato in:", normalizePath(out_dir), "\n\n")