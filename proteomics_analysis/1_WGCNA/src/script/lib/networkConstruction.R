networkConstruction <- function(){
  
  ########################################
  # input parameters
  
  corType <- input_parameter$corType
  power <- input_parameter$power
  networkType <- input_parameter$networkType
  TOMType <- input_parameter$TOMType
  minModuleSize <- input_parameter$minModuleSize
  
  filename_dataInput <- output_file$filename_dataInput
  
  filename_networkConstruction <- output_file$filename_networkConstruction 
  
  filename_moduleSize_barplot <- output_file$filename_moduleSize_barplot
  ########################################
  
  lnames = load(file = filename_dataInput)
  ########################################
  
  powers = c(1:10, seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.80,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  ########################################
  
  net = blockwiseModules(datExpr, 
                         corType = corType,
                         power = power,
                         networkType = networkType,
                         TOMType = TOMType, 
                         minModuleSize = minModuleSize,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = F,
                         verbose = 3)
  ########################################
  
  mergedColors = labels2colors(net$colors)
  
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  ########################################
  
  t <- data.frame(table(moduleColors))
  t <- t[order(t$Freq, decreasing = T),]
  
  pdf(filename_moduleSize_barplot,width=6, height=5)
  par(mar=c(10,5,3,1))
  lab <- paste0(as.character(t$moduleColors),"(",as.character(t$Freq),")")
  barplot(t$Freq, col = as.character(t$moduleColors), 
          ylab = "# Nodes", names.arg = lab, las = 2, cex.names = 1.1)
  dev.off()
  ########################################
  
  save(MEs, moduleLabels, moduleColors, geneTree, file = filename_networkConstruction)

}



