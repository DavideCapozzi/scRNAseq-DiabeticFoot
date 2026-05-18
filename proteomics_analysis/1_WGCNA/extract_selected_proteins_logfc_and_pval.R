library(dplyr)
library(xlsx)

degs_v1 <- read.delim("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V1/Results_V1/DEGs/DEG.txt")
degs_v2 <- read.delim("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V2/Results_V2/DEGs/DEG.txt")
degs_v3 <- read.delim("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/DEGs/DEG.txt")
 
degs_v1_ordered <- degs_v1 %>%
  arrange(FDR)

degs_v2_ordered <- degs_v2 %>%
  arrange(FDR)

degs_v3_ordered <- degs_v3 %>%
  arrange(FDR)

rm(degs_v1, degs_v2, degs_v3)

V1_genes <- c("CCL5", "MMP-1", "CXCL5", "CXCL1", "FCN2", "MCP-4", "C2", "CCL4")
V2_genes <- c("SELL", "CCL18", "CHL1", "F7")
V3_genes <- c("CCL18", "COMP", "FCN2")

V1_ind <- which(degs_v1_ordered$gene %in% V1_genes)
V1_final <- degs_v1_ordered[V1_ind,]

V2_ind <- which(degs_v2_ordered$gene %in% V2_genes)
V2_final <- degs_v2_ordered[V2_ind,]

V3_ind <- which(degs_v3_ordered$gene %in% V3_genes)
V3_final <- degs_v3_ordered[V3_ind,]

write.xlsx(V1_final, "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V1/Results_V1/DEGs/V1_selected_proteins.xlsx", row.names = F)
write.xlsx(V2_final, "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V2/Results_V2/DEGs/V2_selected_proteins.xlsx", row.names = F)
write.xlsx(V3_final, "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/DEGs/V3_selected_proteins.xlsx", row.names = F)

