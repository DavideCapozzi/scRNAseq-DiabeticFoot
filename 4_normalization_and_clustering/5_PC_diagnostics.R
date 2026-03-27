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
N_PCS_USED  <- 6                  # quello attuale
N_PCS_TOTAL <- 50                 # npcs usati in RunPCA
CUTOFFS     <- c(6, 10, 20)       # cutoff da confrontare nei plot [MODIFICATO: rimosso 15]
RES_VALUES  <- c(0.4, 0.5, 0.7)  # resolution per i plot UMAP [AGGIUNTO]


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

# Palette e linetype per 3 cutoff: 6, 10, 20 [MODIFICATO: da 4 a 3 elementi]
cutoff_colors    <- c("firebrick", "#e07b00", "#1a7abd")
cutoff_linetypes <- c("dashed", "dotted", "dotdash")

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
  # marcatori per i cutoff candidati [MODIFICATO: 3 cutoff, 3 colori/linetype]
  geom_vline(xintercept = CUTOFFS + 0.5,
             linetype   = cutoff_linetypes,
             color      = cutoff_colors,
             linewidth  = 0.8) +
  annotate("text",
           x     = CUTOFFS + 0.6,
           y     = max(var_exp) * 0.85,
           label = paste0("PC", CUTOFFS,
                          "\n", sprintf("%.1f%%", cum_var[CUTOFFS])),
           hjust = 0, size = 3,
           color = cutoff_colors) +
  scale_fill_manual(values = c("TRUE" = "#2171b5", "FALSE" = "#bdd7e7")) +
  scale_x_continuous(breaks = c(1, CUTOFFS, seq(10, N_PCS_TOTAL, 10)) |> unique() |> sort()) +
  scale_y_continuous(
    name     = "% Variance explained",
    expand   = expansion(mult = c(0, 0.15)),
    sec.axis = sec_axis(~ . * 4, name = "Cumulative variance (%)")
  ) +
  labs(
    title    = "PCA — Variance explained per PC",
    # Subtitle aggiornato: solo PC6, PC10, PC20 [MODIFICATO]
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
#  3. HEATMAP TOP FEATURE LOADINGS — PC1:20
#     Mostra se le PC tardive (7-20) hanno ancora struttura biologica
#     o sono rumore
#     [MODIFICATO: esteso da PC1:15 a PC1:20 per coerenza con i nuovi cutoff]
# ============================================================
message("[3/6] Heatmap feature loadings PC1-20...")

# DimHeatmap è la funzione nativa Seurat per questo scopo
# Salviamo come PNG aprendo device manualmente (DimHeatmap non è ggplot)
png(file.path(out_dir, "03_dimheatmap_PC1-20.png"),  # [MODIFICATO: filename]
    width = 2400, height = 4000, res = 200)            # [MODIFICATO: h aumentata per 20 PC]
DimHeatmap(seurat_obj,
           dims       = 1:20,   # [MODIFICATO: da 1:15 a 1:20]
           cells      = 500,    # subsample per velocità
           balanced   = TRUE,
           fast       = FALSE)
dev.off()
message("  → salvato: 03_dimheatmap_PC1-20.png")


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
#  5. UMAP SENSITIVITY — confronto UMAP con dims 1:6 / 1:10 / 1:20
#     per ciascuna resolution in c(0.4, 0.5, 0.7)
#     colorato per seurat_clusters attuali
#     NOTA: ricalcola FindNeighbors + RunUMAP in locale,
#           NON sovrascrive l'oggetto originale
#
#  [MODIFICATO]:
#    - PC confrontati: c(6, 10, 20)  (rimosso 15)
#    - Aggiunto loop su resolution: c(0.4, 0.5, 0.7)
#    - Per ogni resolution: griglia 3 UMAP (una per PC cutoff)
#    - Output: un file PNG per resolution → 3 file totali
# ============================================================
message("[5/6] UMAP sensitivity analysis (6 / 10 / 20 PCs  ×  res 0.4 / 0.5 / 0.7)...")

pc_dims <- c(6, 10, 20)  # [MODIFICATO: da c(6, 10, 15) a c(6, 10, 20)]

for (res in RES_VALUES) {
  
  message(sprintf("  Resolution = %.1f", res))
  umap_plots <- list()
  
  for (nd in pc_dims) {
    tmp <- FindNeighbors(seurat_obj, dims = 1:nd, verbose = FALSE)
    tmp <- FindClusters(tmp, resolution = res, verbose = FALSE,
                        graph.name = "RNA_snn")
    tmp <- RunUMAP(tmp, dims = 1:nd, verbose = FALSE,
                   reduction.name = "umap_tmp",
                   reduction.key  = "umaptmp_")
    
    df_umap           <- as.data.frame(tmp@reductions$umap_tmp@cell.embeddings)
    df_umap$cluster   <- tmp$seurat_clusters   # cluster calcolati a questa res
    df_umap$cluster_0 <- seurat_obj$seurat_clusters  # cluster originali (ref)
    
    umap_plots[[as.character(nd)]] <- ggplot(df_umap,
                                             aes(x = umaptmp_1, y = umaptmp_2,
                                                 color = cluster)) +
      geom_point(size = 0.05, alpha = 0.3) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title    = sprintf("UMAP  —  dims = 1:%d  |  res = %.1f", nd, res),
           subtitle = sprintf("Cumulative variance = %.1f%%", cum_var[nd]),
           x = "UMAP 1", y = "UMAP 2") +
      theme_classic(base_size = 11) +
      theme(legend.text  = element_text(size = 7),
            legend.title = element_blank())
    
    rm(tmp); gc()
  }
  
  p_umap_sensitivity <- wrap_plots(umap_plots, ncol = 3) +
    plot_annotation(
      title    = sprintf("UMAP Sensitivity Analysis  —  Resolution = %.1f", res),
      subtitle = sprintf(
        "Clusters calcolati a res=%.1f su UMAP con 1:%d, 1:%d, 1:%d PC rispettivamente",
        res, pc_dims[1], pc_dims[2], pc_dims[3]
      ),
      theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 11))
    )
  
  # Filename con resolution codificata: es. 05_umap_sensitivity_res0.4.png
  fname <- sprintf("05_umap_sensitivity_res%.1f.png", res)
  save_plot(p_umap_sensitivity, fname, w = 18, h = 6)
}


# ============================================================
#  6. ARI/NMI CLUSTER STABILITY — label-permutation-invariant
#     Computes ARI and NMI between clustering at PC6 (reference)
#     and clusterings at PC10, PC20.  [MODIFICATO: rimosso PC15]
#     Per ciascuna resolution in c(0.4, 0.5, 0.7). [AGGIUNTO]
#     ARI = 1 → identical structure; ARI ~ 0 → random agreement
#     Unlike raw concordance, ARI is invariant to cluster renumbering.
# ============================================================
message("[6/6] ARI/NMI cluster stability across PC cutoffs e resolution...")

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

# PC targets: [MODIFICATO da c(10, 15, 20) a c(10, 20)]
pc_targets <- c(10, 20)

# Riferimento: clustering originale dell'oggetto
ref_labels <- as.character(seurat_obj$seurat_clusters)

# Raccoglie risultati per tutte le combinazioni PC × resolution
all_stability <- list()

for (res in RES_VALUES) {
  message(sprintf("  Resolution = %.1f", res))
  
  for (nd in pc_targets) {
    message(sprintf("    Computing ARI/NMI for dims = 1:%d vs 1:%d, res = %.1f ...",
                    N_PCS_USED, nd, res))
    tmp <- FindNeighbors(seurat_obj, dims = 1:nd, verbose = FALSE)
    tmp <- FindClusters(tmp, resolution = res, verbose = FALSE,
                        cluster.name = "clusters_tmp",
                        graph.name   = "RNA_snn")
    cmp_labels <- as.character(tmp$clusters_tmp)
    
    ari_val <- ari_manual(ref_labels, cmp_labels)
    nmi_val <- nmi_manual(ref_labels, cmp_labels)
    
    # contingency matrix (row-normalised) per heatmap tile
    cont <- data.frame(ref = ref_labels, cmp = cmp_labels) |>
      count(ref, cmp) |>
      group_by(ref) |>
      mutate(pct = n / sum(n) * 100) |>
      ungroup() |>
      mutate(
        pc_comparison = sprintf("PC%d vs PC%d", N_PCS_USED, nd),
        resolution    = res,
        ari           = ari_val,
        nmi           = nmi_val
      )
    
    key <- sprintf("nd%d_res%.1f", nd, res)
    all_stability[[key]] <- cont
    rm(tmp); gc()
  }
}

df_stability <- bind_rows(all_stability)

# ── 6a: ARI / NMI summary bar plot (facettato per resolution) ──
df_metrics <- df_stability |>
  distinct(pc_comparison, resolution, ari, nmi) |>
  tidyr::pivot_longer(cols = c(ari, nmi),
                      names_to  = "metric",
                      values_to = "value") |>
  mutate(metric     = toupper(metric),
         label      = sprintf("%.3f", value),
         resolution = factor(sprintf("res = %.1f", resolution),
                             levels = sprintf("res = %.1f", RES_VALUES)))

p_ari <- ggplot(df_metrics, aes(x = pc_comparison, y = value, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = label),
            position = position_dodge(0.7),
            vjust = -0.4, size = 3.2) +
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
  facet_wrap(~ resolution, ncol = 3) +   # [AGGIUNTO: una colonna per res]
  labs(
    title    = "Cluster stability: ARI and NMI (label-permutation invariant)",
    subtitle = sprintf("Reference: clustering at dims = 1:%d  |  Resolutions = %s",
                       N_PCS_USED, paste(RES_VALUES, collapse = " / ")),
    x = NULL, y = "Score (0 = random, 1 = identical)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold"))

save_plot(p_ari, "06a_ARI_NMI_summary.png", w = 14, h = 6)  # [MODIFICATO: w allargato]

# ── 6b: contingency heatmaps — una griglia per resolution ───
# Per ogni resolution: heatmap per ogni PC target (PC10, PC20)
# [MODIFICATO: rimosso PC15, aggiunto loop su resolution]

for (res in RES_VALUES) {
  
  heatmap_list <- lapply(pc_targets, function(nd) {
    sub <- df_stability |>
      filter(pc_comparison == sprintf("PC%d vs PC%d", N_PCS_USED, nd),
             abs(resolution - res) < 1e-9)
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
        x = sprintf("Clusters (dims = 1:%d, res = %.1f)", nd, res),
        y = sprintf("Clusters (dims = 1:%d, ref)", N_PCS_USED)
      ) +
      theme_classic(base_size = 10) +
      theme(axis.text       = element_text(size = 7),
            legend.position = "right")
  })
  
  p_heatmaps <- wrap_plots(heatmap_list, ncol = 2) +  # [MODIFICATO: ncol=2 per 2 PC target]
    plot_annotation(
      title = sprintf(
        "Contingency heatmaps: reference clustering (PC%d) vs PC cutoffs  |  res = %.1f",
        N_PCS_USED, res
      ),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  fname <- sprintf("06b_contingency_heatmaps_res%.1f.png", res)
  save_plot(p_heatmaps, fname, w = 16, h = 8)
}


# ============================================================
#  SUMMARY REPORT
#  [MODIFICATO]: rimosso PC15, aggiunto loop su resolution
# ============================================================

# Estrai ARI per ogni combinazione PC × resolution
get_ari <- function(pc_nd, res) {
  df_stability |>
    filter(pc_comparison == sprintf("PC%d vs PC%d", N_PCS_USED, pc_nd),
           abs(resolution - res) < 1e-9) |>
    pull(ari) |>
    unique()
}

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
cat(sprintf("  Variance @ PC20 : %.2f%%\n", cum_var[20]))
cat(sprintf("  Variance @ PC50 : %.2f%%\n", cum_var[50]))
cat("----------------------------------------------\n")
cat("  Cluster stability (ARI, label-invariant):\n")
for (res in RES_VALUES) {
  cat(sprintf("  Resolution = %.1f:\n", res))
  for (nd in pc_targets) {
    ari_val <- get_ari(nd, res)
    cat(sprintf("    PC%d vs PC%d : ARI = %.3f\n", N_PCS_USED, nd, ari_val))
  }
}
cat("  Interpretation: ARI > 0.8 = strong stability\n")
cat("                  ARI 0.6-0.8 = moderate\n")
cat("                  ARI < 0.6 = substantial restructuring\n")
cat("══════════════════════════════════════════════\n")
cat("  Output saved in:", normalizePath(out_dir), "\n\n")