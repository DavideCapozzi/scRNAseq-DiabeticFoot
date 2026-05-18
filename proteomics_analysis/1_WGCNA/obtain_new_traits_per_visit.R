library(xlsx)

new_clinical_data <- read.xlsx("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/Clinical_Data_patients_new.xlsx", sheetIndex = 1)
colnames(new_clinical_data) <- gsub("\\.", " ", colnames(new_clinical_data))
#rownames(new_clinical_data) <- new_clinical_data$`Proteomics ID `

V1 <- read.table("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V1/matrix/Traits.txt", header = T, sep = "\t", quote = "", check.names = F)
V1_samples <- V1$ProteiomicsID
V1_new_index <- which(new_clinical_data$`Proteomics ID ` %in% V1$ProteiomicsID)
V1_new_clinical_data <- new_clinical_data[V1_new_index,]

V2 <- read.table("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V2/matrix/Traits.txt", header = T, sep = "\t", quote = "", check.names = F)
V2_samples <- V2$ProteiomicsID
V2_new_index <- which(new_clinical_data$`Proteomics ID ` %in% V2$ProteiomicsID)
V2_new_clinical_data <- new_clinical_data[V2_new_index,]
