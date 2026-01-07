# =========================================================
# LIBRERIE
# =========================================================
library(Seurat)
library(dplyr)
library(openxlsx)
library(ggplot2)

set.seed(123)

# =========================================================
# CONFIGURAZIONE
# =========================================================
output_dir <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/5_annotation/emelie_markers_featureplots/"
dir.create(output_dir)
seurat_obj <- readRDS("/home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/res_clustering/seurat_res_0.7.rds")
resolution <- 0.7
genes_excel <- "/home/fdann/Desktop/proj/sc-res/2_new_analysis/5_annotation/emelie_markers.xlsx"
heatmap_colors <- c("#336699", "white", "#CC0000")

# =========================================================
# IDENTIFICAZIONE DEI CLUSTER
# =========================================================
Idents(seurat_obj) <- paste0("RNA_snn_res.", resolution)

cat("\n==============================\n")
cat("HEATMAP DEI MARCATORI PER CLUSTER - RESOLUTION =", resolution, "\n")
cat("==============================\n")

# =========================================================
# CARICAMENTO LISTA DI GENI DA EXCEL
# =========================================================
cat("\nCaricamento lista di geni da Excel...\n")

# Il file deve contenere una singola colonna con i nomi dei geni
genes_df <- read.xlsx(genes_excel, sheet = 2, colNames = TRUE)

# Cerca una colonna chiamata "gene" o usa la prima colonna
if ("gene" %in% tolower(names(genes_df))) {
  colname <- names(genes_df)[tolower(names(genes_df)) == "gene"]
  custom_genes <- genes_df[[colname]]
} else {
  custom_genes <- genes_df[[1]]
}

custom_genes <- unique(custom_genes)
cat("Numero totale di geni caricati:", length(custom_genes), "\n")

# Filtra solo i geni presenti nel Seurat object
present_genes <- custom_genes[custom_genes %in% rownames(seurat_obj)]
missing_genes <- setdiff(custom_genes, present_genes)

cat("Geni trovati nel Seurat object:", length(present_genes), "\n")
cat("Geni mancanti:", length(missing_genes), "\n")

if (length(missing_genes) > 0) {
  cat("⚠️  ATTENZIONE: alcuni geni non trovati nel dataset:\n")
  print(missing_genes)
}

if (length(present_genes) == 0) {
  stop("❌ Nessuno dei geni caricati è presente nel Seurat object.")
}

# =========================================================
# SCALING DEI GENI PER LA HEATMAP
# =========================================================
cat("\nScalatura dei geni selezionati per la heatmap...\n")

tmp_obj <- ScaleData(seurat_obj, features = present_genes, verbose = FALSE)

# =========================================================
# CREAZIONE HEATMAP
# =========================================================
cat("Plot della heatmap...\n")

heatmap_plot <- DoHeatmap(
  tmp_obj,
  features = present_genes,
  group.by = paste0("RNA_snn_res.", resolution),
  size = 4
) + 
  scale_fill_gradientn(colors = heatmap_colors) +
  ggtitle("Literature markers heatmap") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Mostra la heatmap
print(heatmap_plot)

# =========================================================
# SALVATAGGIO DELLA HEATMAP
# =========================================================
output_file <- file.path(output_dir, "heatmap_geni_excel_markers.png")
ggsave(output_file, plot = heatmap_plot, width = 12, height = 10)
cat("\n✅ Heatmap salvata in:", output_file, "\n")

# =========================================================
# CREAZIONE FEATUREPLOT PER OGNI GENE
# =========================================================
cat("\nCreazione dei FeaturePlot per ogni gene...\n")

# Directory dove salvare i plot dei FeaturePlot
featureplot_dir <- file.path(output_dir, "featureplots")
if (!dir.exists(featureplot_dir)) {
  dir.create(featureplot_dir, recursive = TRUE)
}

# Loop su ogni gene presente
for (gene in present_genes) {
  cat("Plotting gene:", gene, "...\n")
  
  p <- FeaturePlot(
    seurat_obj,
    features = gene,
    reduction = "umap",
    cols = c("lightgrey", "#CC0000")
  ) +
    ggtitle(paste("FeaturePlot:", gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Salvataggio del plot
  output_fp <- file.path(featureplot_dir, paste0("FeaturePlot_", gene, ".png"))
  ggsave(output_fp, plot = p, width = 6, height = 5)
}

cat("\n✅ Tutti i FeaturePlot salvati in:", featureplot_dir, "\n")

