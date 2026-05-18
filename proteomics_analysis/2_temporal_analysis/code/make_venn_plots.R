library(ggvenn)

venn_list <- list(
  V1 = c("CCL5", "MMP-1", "CXCL5", "CXCL1", "FCN2", "MCP-4", "C2", "CCL4"),
  V2 = c("SELL", "CCL18", "CHL1", "F7"),
  V3 = c("CCL18", "COMP", "FCN2"))

  p_venn <- ggvenn(
    venn_list, 
    fill_color = c("green", "blue2", "turquoise"),
    stroke_size = 0.5, 
    show_percentage = F ,
    set_name_size = 4
  ) + ggtitle("Proteins overlap between visits")+
    theme(plot.title = element_text(hjust = 0.5, size = 14)) 
  
pdf("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/temporal_analysis/venn.pdf", width = 7, height = 5)
print(p_venn)
dev.off()
  

 


common_all <- intersect(intersect(venn_list$V1, venn_list$V2), venn_list$V3)

v1_v2_only <- setdiff(intersect(venn_list$V1, venn_list$V2), venn_list$V3)
v1_v3_only <- setdiff(intersect(venn_list$V1, venn_list$V3), venn_list$V2)
v2_v3_only <- setdiff(intersect(venn_list$V2, venn_list$V3), venn_list$V1)
  

v1_only <- setdiff(venn_list$V1, union(venn_list$V2, venn_list$V3))
v2_only <- setdiff(venn_list$V2, union(venn_list$V1, venn_list$V3))
v3_only <- setdiff(venn_list$V3, union(venn_list$V1, venn_list$V2))
  
cat("--- PROTEINE ESCLUSIVE ---\n")
cat("Solo V1:    ", paste(v1_only, collapse = ", "), "\n")
cat("Solo V2:    ", paste(v2_only, collapse = ", "), "\n")
cat("Solo V3:    ", paste(v3_only, collapse = ", "), "\n\n")

cat("--- INTERSEZIONI DOPPIE ---\n")
cat("V1 & V2:    ", paste(v1_v2_only, collapse = ", "), "\n")
cat("V1 & V3:    ", paste(v1_v3_only, collapse = ", "), "\n")
cat("V2 & V3:    ", paste(v2_v3_only, collapse = ", "), "\n\n")

cat("--- INTERSEZIONE TRIPLA ---\n")
cat("Tutte (V1, V2, V3): ", paste(common_all, collapse = ", "), "\n")
