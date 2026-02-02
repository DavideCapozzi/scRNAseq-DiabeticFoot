library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Imposta la working directory
setwd("/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/seurat_res_0.7")

# Carica il Seurat object
seu <- readRDS("seurat_res_0.7.rds")

# Trova i marker per cluster
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers <- markers %>% arrange(cluster, desc(avg_log2FC))

# Salva la tabella completa
write.csv(markers, "/home/fdann/Desktop/proj/sc-res/2_new_analysis/5_annotation/upregulated_markers/markers_per_cluster_res0.7.csv", row.names = FALSE)

# Numero di geni upregolati da includere nella heatmap
top_n <- 3

# Seleziona i top marker per ogni cluster
top_markers <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = top_n, wt = avg_log2FC)

# Ordina i geni per cluster
top_markers <- top_markers %>% arrange(cluster, desc(avg_log2FC))

genes_for_heatmap <- unique(top_markers$gene)

# Heatmap: righe = geni up, colonne = cluster
pdf("/home/fdann/Desktop/proj/sc-res/2_new_analysis/5_annotation/upregulated_markers/heatmap_top_markers_clusterwise_res0.7.pdf", width = 14, height = 10)
DoHeatmap(
  object = seu,
  features = genes_for_heatmap,
  group.by = "seurat_clusters",
  size = 3
) + 
  ggtitle("Top upregulated genes per cluster (resolution 0.7)")
dev.off()

cat("Analisi completata: marker salvati e heatmap generata.\n")

