rm(list = ls())
options(stringAsFactors = F)

###########################################
library(stringr)
library(pheatmap)
library(xlsx)
##########################################
path <- "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/"
visit <- "DiabeticFoot_V1"
##########################################
dirRes <- paste0(path, visit, "/Results_V1/")

if(!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory", dirRes,"already exists"))
}

DEGs_folder <- "DEGs"

dirDEGs <- paste0(dirRes,DEGs_folder,"/")

if(!dir.exists(dirDEGs)){
  dir.create(dirDEGs)
}else{
  print(paste("The directory", dirRes,"already exists"))
}

########################################## 
filename_in <- paste0(path, visit, "/matrix", "/matrix.txt")
condition_table <- read.xlsx("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/condition_sample_match.xlsx", sheetIndex = 1 )

#output
file_DEG <- paste0(dirDEGs, "DEG.txt")
filename_matrix_DEG <- paste0(dirDEGs, "matrix_DEG.txt")
filename_heatmap <- paste0(dirDEGs, "heatmap.pdf")
###########################################
#parameters setting

prc <- 0
soglia_fc <- 1
soglia_pval <- 1
##########################################
#STEP 1: import data
tmp<- read.table(filename_in, header = T, check.names = F, row.names = 1, sep = "\t")
#classes <- sapply(tmp, class)
#tmp<- read.table(filename_in, header = T, check.names = F, row.names = 1, sep = "\t", quote = "", colClasses = classes)

geni <- rownames(tmp)
pz <- colnames(tmp)

not_healed <- condition_table[grep("not_healed", condition_table$condition), ]
healed <- condition_table[!grepl("not_healed", condition_table$condition), ]

not_healed_present <- which(colnames(tmp) %in% not_healed$Proteomics.ID)
healed_present <- which(colnames(tmp) %in%  healed$Proteomics.ID )

dataN <- tmp[,not_healed_present]
dataC <- tmp[,healed_present]
data <- tmp
# exp_N <- "TCGA-\\w+-\\w+-1\\d"
# exp_C <- "TCGA-\\w+-\\w+-0\\d"
# 
# head(pz)
# 
# pzN <- grep(exp_N, pz, value = T)
# pzC <- grep(exp_C, pz, value = T)

# name_N <- str_extract(pzN, "TCGA-\\w+-\\w+")
# name_C <- str_extract(pzC, "TCGA-\\w+-\\w+")
# 
# pzC <- pzC[!duplicated(name_C)]
# pzN <- pzN[!duplicated(name_N)]
# 
# common_names <- intersect(name_N, name_C)

#pzC <- 1:length(common_names)
#pzN <- 1:length(common_names)

#for (k in 1:length(common_names)) {
  #pzC[k] <- grep(name_C[k], pzC, value = T)
  #pzN[k] <- grep(name_N[k], pzN, value = T)
#}

# pzN_com <- unlist(lapply(common_names, function(x){grep(x, pzN, value=TRUE)}))
# pzC_com <- unlist(lapply(common_names, function(x){grep(x, pzC, value=TRUE)}))

# dataN <- tmp[, pzN_com]
# dataC <- tmp[, pzC_com]
#data <- cbind(dataN, dataC)

#pz_com <- colnames(data)


# dataN <- log2(dataN + 1)
# dataC <- log2(dataC + 1)
# data <- log2(data + 1)

#rm(tmp,pz,name_N, name_C, pzN, pzC, common_names, classes)
#pre processing #calculate iqr for all the matrix
variation <- apply(data, 1, IQR, type=5)
soglia_prc <- quantile(variation, prc)
ind <- which(variation <= soglia_prc)
dataN <- dataN[-ind,]
dataC <- dataC[-ind,]
data <- data[-ind,]
geni <- geni[-ind]
rm(ind)

hist(variation,
      main ="IQR freq dist", breaks=100,
      xlab="IQR", ylab="Genes freq", col="green") 

abline(v = soglia_prc, lty=2, lwd=4, col="red")
#filtering
logFC <- rowMeans(dataC) - rowMeans(dataN)

hist(logFC,
     main ="logFC freq distribution", breaks=100,
     xlab="logFC", ylab="Genes freq", col="red") 
abline(v = c(-log2(soglia_fc), log2(soglia_fc)), lty=2, lwd=4, col="grey")

ind <- which(abs(logFC) < log2(soglia_fc))
if(length(ind)>0) {
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  geni <- geni[-ind]
  logFC <- logFC[-ind]
}

rm(ind)  
#pvalue
N <- ncol(dataN)
M <- ncol(dataC)

pval <- apply(data, 1, function(x) {
  res <- t.test(x[1:N], x[(N+1):(M+N)], paired = F)
  pval <- res$p.value
  return(pval)
})
#corrected pvalue
pval_adj <- p.adjust(pval, method = "fdr")

#filtro su pvalue
ind <- which(pval > soglia_pval)

if(length(ind)>0) {
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  geni <- geni[-ind]
  logFC <- logFC[-ind]
  pval <- pval[-ind]
}

rm(ind)
#plots
#volcano plot
plot(logFC, -log10(pval),
     main = "volcano plot",
     xlim = c(-8,8),
     ylim = c(0,70),
     xlab ="log2FC",
     ylab="-log10 pval")

abline(h = -log10(soglia_pval), lty=2, lwd=4, col="red")
abline(v = c(-log2(soglia_fc), log2(soglia_fc)), lty=2, lwd=4, col="green")

#box plot
#ordered_FC <- logFC[order(logFC, decreasing = T)]
gene_id <- "CCR1"
ind <- which(geni %in% gene_id) 
df <- data.frame(normal = t(dataN[ind,]),
                  cancer = t(dataC[ind,]),
                  row.names = NULL)
colnames(df) <- c("not_healed","healed")

pdf("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/res_validation/interactors_analysis/CCL3_boxplot_notch.pdf")
boxplot(df,
        main = paste0(gene_id, ",","p-value= ",
                      format(pval[ind], digits =2)),
        ylab = "gene expression value",
        xlab= "condition",
        col = c("green", "orange"),
        notch = T, #confidence interval for 0.05 (i can say the gene is significant with a confidence of at least 0.05)
        pars = list(boxwex = 0.3, staplewex = 0.6))
dev.off()
rm(ind)
#exporting results
direction <- ifelse(logFC>0, "up", "down")
results <- data.frame(geni, pval, pval_adj, logFC, direction)
barplot(table(direction))
write.table(results, file=file_DEG, row.names=F,
            col.names=c("gene", "pval","FDR", "log2-FC", "direction"),
            sep="\t", quote=FALSE)
#heatmap
test <- grepl('TCGA-\\w+-\\w+-1\\d', colnames(data))
samples <- ifelse(test, "normal", "cancer")

annotation <- data.frame(SampleType=samples)
rownames(annotation) <- colnames(data)

vect_color <- c("green", "orange")
names(vect_color) <- unique(samples)

annotation_colors <- list(SampleType = vect_color)

pheatmap(data, scale ="row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         color = colorRampPalette(colors=c("blue", "blue3","black","yellow3", "yellow")) (100),
         show_rownames = F,
         show_colnames = F,
         cutree_cols =2,
         cutree_rows =2,
         width=10, height=10,
         file = filename_heatmap
)

