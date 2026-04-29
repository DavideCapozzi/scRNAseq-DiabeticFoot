library(xlsx)

ccl4_interactors <- read.xlsx("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/interactors_analysis/CCL4_interactors.xlsx", sheetIndex = 1)
cor_table <- read.table("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/analysis_joined_WGCNA_paola/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V1/Results_V1/txtFile/cytoscapeInput_nodes.txt", sep = "\t", quote = "", check.names = F, header = T)

x <- which(cor_table$nodeName %in% ccl4_interactors$Interactor)
interactors <- cor_table[x,] # just ccl3 found in proteomics (module turquoise)
##################

degs_sc <- read.xlsx("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusts.xlsx", sheetIndex =1)

x <- which(degs_sc$Gene %in% ccl4_interactors$Interactor)
interactors_sc <- degs_sc[x,] 

ind <-which(interactors_sc$Gene == "CCL4")
interactors_sc <- interactors_sc[-ind,]

unique(interactors_sc$Gene) # 10 interactors founds in degs 

write.xlsx(interactors_sc, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/interactors_analysis/CCL4_interactors_in_sc.xlsx", 
           row.names = F)

