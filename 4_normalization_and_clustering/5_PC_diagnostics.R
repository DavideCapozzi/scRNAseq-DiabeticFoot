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

# ── output dir ───────────────────────────────────────────────
out_dir <- "./PC_diagnostics"
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
#  6. ARI/NMI CLUSTER STABILITY — label-permutation-invariant
#     Computes ARI and NMI between clustering at PC6 (reference)
#     and clusterings at PC10, PC15, PC20.
#     ARI = 1 → identical structure; ARI ~ 0 → random agreement
#     Unlike raw concordance, ARI is invariant to cluster renumbering.
# ============================================================
message("[6/6] ARI/NMI cluster stability across PC cutoffs...")

# bluster::pairwiseRand requires Bioconductor; use fossil::adj.rand.index
# as lightweight alternative, falling back to manual computation if needed.
ari_manual <- function(a, b) {
  # Computes ARI from two factor/character vectors via contingency table.
  # Hubert & Arabie (1985) formulation.
  a <- as.integer(factor(a))
  b <- as.integer(factor(b))
  tab <- table(a, b)
  n   <- sum(tab)
  sa  <- sum(choose(rowSums(tab), 2))
  sb  <- sum(choose(colSums(tab), 2))
  sc  <- sum(choose(tab, 2))
  expected <- sa * sb / choose(n, 2)
  denom    <- 0.5 * (sa + sb) - expected
  if (denom == 0) return(1.0)
  (sc - expected) / denom
}

nmi_manual <- function(a, b) {
  # Normalized Mutual Information (arithmetic mean normalisation).
  a  <- as.integer(factor(a))
  b  <- as.integer(factor(b))
  n  <- length(a)
  tab <- table(a, b) / n
  pa  <- rowSums(tab)
  pb  <- colSums(tab)
  # mutual information
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      pij <- tab[i, j]
      if (pij > 0) mi <- mi + pij * log(pij / (pa[i] * pb[j]))
    }
  }
  ha <- -sum(pa[pa > 0] * log(pa[pa > 0]))
  hb <- -sum(pb[pb > 0] * log(pb[pb > 0]))
  denom <- 0.5 * (ha + hb)
  if (denom == 0) return(1.0)
  mi / denom
}

ref_labels <- as.character(seurat_obj$seurat_clusters)
pc_targets  <- c(10, 15, 20)
stability_results <- list()

for (nd in pc_targets) {
  message(sprintf("  Computing ARI/NMI for dims = 1:%d vs 1:%d ...", N_PCS_USED, nd))
  tmp <- FindNeighbors(seurat_obj, dims = 1:nd, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.7, verbose = FALSE,
                      cluster.name = "clusters_tmp",
                      graph.name   = "RNA_snn")
  cmp_labels <- as.character(tmp$clusters_tmp)
  
  ari_val <- ari_manual(ref_labels, cmp_labels)
  nmi_val <- nmi_manual(ref_labels, cmp_labels)
  
  # contingency matrix (row-normalised) for heatmap tile
  cont <- data.frame(ref = ref_labels, cmp = cmp_labels) |>
    count(ref, cmp) |>
    group_by(ref) |>
    mutate(pct = n / sum(n) * 100) |>
    ungroup() |>
    mutate(pc_comparison = sprintf("PC%d vs PC%d", N_PCS_USED, nd),
           ari = ari_val, nmi = nmi_val)
  
  stability_results[[as.character(nd)]] <- cont
  rm(tmp); gc()
}

df_stability <- bind_rows(stability_results)

# ── 6a: ARI / NMI summary bar plot ──────────────────────────
df_metrics <- df_stability |>
  distinct(pc_comparison, ari, nmi) |>
  tidyr::pivot_longer(cols = c(ari, nmi),
                      names_to  = "metric",
                      values_to = "value") |>
  mutate(metric = toupper(metric),
         label  = sprintf("%.3f", value))

p_ari <- ggplot(df_metrics, aes(x = pc_comparison, y = value, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = label),
            position = position_dodge(0.7),
            vjust = -0.4, size = 3.5) +
  geom_hline(yintercept = c(0.6, 0.8),
             linetype = c("dashed", "dotted"),
             color    = c("#e07b00", "#2171b5"),
             linewidth = 0.7) +
  annotate("text", x = 0.55, y = 0.62,
           label = "ARI = 0.6\n(moderate)", color = "#e07b00",
           size = 3, hjust = 0) +
  annotate("text", x = 0.55, y = 0.82,
           label = "ARI = 0.8\n(strong)", color = "#2171b5",
           size = 3, hjust = 0) +
  scale_fill_manual(values = c("ARI" = "#2171b5", "NMI" = "#cb181d"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Cluster stability: ARI and NMI (label-permutation invariant)",
    subtitle = sprintf("Reference: clustering at dims = 1:%d  |  Resolution = 0.7", N_PCS_USED),
    x = NULL, y = "Score (0 = random, 1 = identical)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

save_plot(p_ari, "06a_ARI_NMI_summary.png", w = 9, h = 6)

# ── 6b: contingency heatmaps per ogni confronto ─────────────
heatmap_list <- lapply(pc_targets, function(nd) {
  sub <- df_stability |>
    filter(pc_comparison == sprintf("PC%d vs PC%d", N_PCS_USED, nd))
  ari_val <- unique(sub$ari)
  nmi_val <- unique(sub$nmi)
  
  ggplot(sub, aes(x = cmp, y = ref, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct > 3, sprintf("%.0f%%", pct), "")),
              size = 2.5, color = "white") +
    scale_fill_gradient(low = "white", high = "#08306b",
                        name = "% cells\n(row)", limits = c(0, 100)) +
    labs(
      title    = sprintf("PC%d vs PC%d", N_PCS_USED, nd),
      subtitle = sprintf("ARI = %.3f  |  NMI = %.3f", ari_val, nmi_val),
      x = sprintf("Clusters (dims = 1:%d)", nd),
      y = sprintf("Clusters (dims = 1:%d)", N_PCS_USED)
    ) +
    theme_classic(base_size = 10) +
    theme(axis.text      = element_text(size = 7),
          legend.position = "right")
})

p_heatmaps <- wrap_plots(heatmap_list, ncol = 3) +
  plot_annotation(
    title = "Contingency heatmaps: reference clustering (6 PC) vs increasing PC cutoffs",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

save_plot(p_heatmaps, "06b_contingency_heatmaps.png", w = 20, h = 8)


# ============================================================
#  SUMMARY REPORT
# ============================================================
ari_10 <- df_stability |> filter(pc_comparison == "PC6 vs PC10") |>
  pull(ari) |> unique()
ari_15 <- df_stability |> filter(pc_comparison == "PC6 vs PC15") |>
  pull(ari) |> unique()
ari_20 <- df_stability |> filter(pc_comparison == "PC6 vs PC20") |>
  pull(ari) |> unique()

cat("\n")
cat("══════════════════════════════════════════════\n")
cat("  PC DIAGNOSTICS — SUMMARY\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Cells            : %d\n", ncol(seurat_obj)))
cat(sprintf("  Genes            : %d\n", nrow(seurat_obj)))
cat(sprintf("  Variable features: %d\n", length(VariableFeatures(seurat_obj))))
cat(sprintf("  PCs computed     : %d\n", N_PCS_TOTAL))
cat("----------------------------------------------\n")
cat(sprintf("  Variance @ PC6  : %.2f%%\n", cum_var[6]))
cat(sprintf("  Variance @ PC10 : %.2f%%\n", cum_var[10]))
cat(sprintf("  Variance @ PC15 : %.2f%%\n", cum_var[15]))
cat(sprintf("  Variance @ PC20 : %.2f%%\n", cum_var[20]))
cat(sprintf("  Variance @ PC50 : %.2f%%\n", cum_var[50]))
cat("----------------------------------------------\n")
cat("  Cluster stability (ARI, label-invariant):\n")
cat(sprintf("    PC6 vs PC10 : ARI = %.3f\n", ari_10))
cat(sprintf("    PC6 vs PC15 : ARI = %.3f\n", ari_15))
cat(sprintf("    PC6 vs PC20 : ARI = %.3f\n", ari_20))
cat("  Interpretation: ARI > 0.8 = strong stability\n")
cat("                  ARI 0.6-0.8 = moderate\n")
cat("                  ARI < 0.6 = substantial restructuring\n")
cat("══════════════════════════════════════════════\n")
cat("  Output saved in:", normalizePath(out_dir), "\n\n")