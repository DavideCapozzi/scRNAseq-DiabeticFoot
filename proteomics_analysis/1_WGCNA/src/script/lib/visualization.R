visualization <- function(){
  
  ########################################
  # input parameters
  
  power <- input_parameter$power
  
  trait_name <- input_parameter$trait_name
  
  filename_dataInput <- output_file$filename_dataInput
  
  filename_networkConstruction <- output_file$filename_networkConstruction
  ########################################
  
  lnames = load(file = filename_dataInput)
  
  lnames = load(file = filename_networkConstruction)
  ######################################## 
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = power);
  # # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  # plotTOM = dissTOM^7;
  # diag(plotTOM) = NA;
  # sizeGrWindow(9,9)
  # TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  ######################################## 
  
  nSelect = 400
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  sizeGrWindow(9,9)
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  ######################################## 
  
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
    trait = as.data.frame(datTraits[, trait_name])
  names(trait) = trait_name

  MET = orderMEs(cbind(MEs, trait))

  par(cex = 0.8)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,4.5), marHeatmap = c(6,5,1,2), cex.lab = 0.9, xLabelsAngle
                        = 90)

  # sizeGrWindow(6,6);
  # par(cex = 1.0)
  # plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
  #                       plotHeatmaps = FALSE)
  # # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  # par(cex = 1.0)
  # plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
  #                       plotDendrograms = FALSE, xLabelsAngle = 90)
  
  
}



