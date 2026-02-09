library(xlsx)

df <- read.xlsx("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusters.xlsx", sheetIndex = 1)

# 2. Definisci la mappatura (basata sul contenuto di celltype_id.txt)
new_names <- c(
  "0" = "T cells (cytotoxic gamma delta)", 
  "1" = "CD4 T cells (CD4 memory/activated Th2 cells)",
  "2" = "CD4 T cells (naive/CM)",
  "3" = "NK cells (CD16)",
  "4" = "CD4 T cells (activated CD4 memory T cells)",
  "5" = "CD8 T cells (EM)",
  "6" = "CD4 T cells (naive/CM)",
  "7" = "T cells (cytotoxic gamma delta)",
  "8" = "CD4 T cells (naive/CM)",
  "9" = "B cells (activated/pre-plasmablasts)", 
  "10" = "Classical monocytes",
  "11" = "B cells (naive/transitional)",
  "12" = "T cells (NKT-like gamma delta)",
  "13" = "CD4/CD8 T cells (memory T cells)",
  "14" = "B cells (naive/transitional)", 
  "15" = "Nonclassical monocytes",
  "16" = "CD4 T cells (naive/CM)",
  "17" = "Intermediate monocytes",
  "18" = "pDCs & cDC2"
)

# 3. Aggiungi la colonna 'cell type'
# Convertiamo la colonna 'cluster' in formato testo per farla corrispondere agli ID del vettore
df$cell_type <- new_names[as.character(df$cluster)]

# 4. (Opzionale) Salva il risultato in un nuovo file CSV
write.xlsx(df, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/scRNAseq-DiabeticFoot/6_differential_expression_analysis/results/DE_with_cell_types.xlsx", row.names = FALSE)
