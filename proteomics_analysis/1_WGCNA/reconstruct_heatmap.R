library(WGCNA)

# ==========================================
# 1. CARICAMENTO DATI
# ==========================================

# Inserisci qui i percorsi ai file generati dalla Parte 1
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/txtFile/"
file_cor  <- paste0(base_dir, "moduleTraitCorV3.txt")
file_pval <- paste0(base_dir, "moduleTraitPvalueV3.txt")

# Leggiamo le tabelle preservando i nomi di riga e colonna
# check.names = FALSE è importante per mantenere caratteri speciali nei nomi dei tratti se presenti
modTraitCor    <- as.matrix(read.table(file_cor, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
modTraitPvalue <- as.matrix(read.table(file_pval, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))

# ==========================================
# 2. PREPARAZIONE DEL PLOT
# ==========================================

# Ricostruiamo il testo che va dentro le celle: "Correlazione \n (P-value)"
# Usiamo le stesse cifre significative dello script originale (2 per cor, 1 per pval)
textMatrix <- paste(signif(modTraitCor, 2), "\n(",
                    signif(modTraitPvalue, 1), ")", sep = "")

# È fondamentale ridare le dimensioni corrette alla matrice di testo
dim(textMatrix) <- dim(modTraitCor)

# ==========================================
# 3. GENERAZIONE HEATMAP (Con Fix Margini)
# ==========================================

#output_pdf <- "G:/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/figure/Heatmap_Ricostruita.pdf"
output_pdf <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/figure/Heatmap_Ricostruita.pdf"

# Aumentiamo slightly width e height per dare respiro
pdf(output_pdf, width = 10, height = 8)

# --- IL FIX PER LE ETICHETTE ---
# par(mar = c(bottom, left, top, right))
# Aumentiamo il primo valore (bottom) da 6 a 12 (o più) per far stare le etichette verticali/ruotate
par(mar = c(12, 9, 3, 3))

labeledHeatmap(Matrix = modTraitCor,
               xLabels = colnames(modTraitCor),   # Nomi dei tratti (Asse X)
               yLabels = rownames(modTraitCor),   # Nomi dei moduli (Asse Y)
               ySymbols = rownames(modTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,           # Il testo ricostruito con i p-values corretti
               setStdMargins = FALSE,             # Disabilitiamo i margini standard per usare i nostri custom
               cex.text = 0.5,                    # Grandezza testo nelle celle
               cex.lab.x = 0.7,                   # Grandezza etichette asse X (riduci se si sovrappongono)
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

message("Heatmap salvata correttamente in: ", getwd(), "/", output_pdf)