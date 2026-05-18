dataInput <- function(){
  
  ########################################
  # input parameters
  
  abline_h <- input_parameter$abline_h
  
  data <- input_file$data
  traits <- input_file$traits
  
  filename_dataInput <- output_file$filename_dataInput
  ########################################
  
  datExpr0 = as.data.frame(t(data))
  gsg = goodSamplesGenes(datExpr0, verbose = 3)
  
  if(!gsg$allOK){
    
    if( sum(!gsg$goodGenes) > 0 ){
      
      printFlush(paste("Removing genes:", 
                       paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
      
    }
    
    if( sum(!gsg$goodSamples) > 0){
      
      printFlush(paste("Removing samples:", 
                       paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
      
    }
    
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
    
  }
  ########################################
  # 
  sampleTree = hclust(dist(datExpr0), method = "average")

  sizeGrWindow(12,9)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree,
       main = "Sample clustering to detect outliers",
       sub ="",
       xlab = "",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 2)

  abline(h = abline_h, col = "red");

  clust = cutreeStatic(sampleTree, cutHeight = abline_h, minSize = 10)

  # keepSamples = (clust == 1)
  # datExpr = datExpr0[keepSamples, ]
  datExpr = datExpr0

  datTraits <- traits[rownames(datExpr),]
  ########################################
  
  sampleTree = hclust(dist(datExpr), method = "average")
  
  traitColors = numbers2colors(datTraits, signed = FALSE)
  
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(datTraits), 
                      dendroLabels = FALSE,
                      main = "Sample dendrogram and trait heatmap")
  ########################################
  
  save(datExpr, datTraits, file = filename_dataInput)

}