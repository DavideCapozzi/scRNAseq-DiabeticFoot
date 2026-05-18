exportNetwork <- function(){
  
  ########################################
  # input parameters
  
  power <- input_parameter$power
  TOMType <- input_parameter$TOMType
  
  module <- input_parameter$module
  
  filename_dataInput <- output_file$filename_dataInput
  
  filename_networkConstruction <- output_file$filename_networkConstruction
  
  filename_cytoscapeInput_edges <-   output_file$filename_cytoscapeInput_edge
  filename_cytoscapeInput_nodes <- output_file$filename_cytoscapeInput_nodes
  ########################################
  
  lnames = load(file = filename_dataInput)
  
  lnames = load(file = filename_networkConstruction)
  ########################################   
  
  TOM = TOMsimilarityFromExpr(datExpr, power = power, TOMType = TOMType); 
  
  #modules = module
  modules = unique(moduleColors)
  
  nodes = names(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modNodes = nodes[inModule];
  
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modNodes, modNodes)
  ######################################## 
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = filename_cytoscapeInput_edges, 
                                 nodeFile = filename_cytoscapeInput_nodes,
                                 weighted = TRUE,
                                 threshold = 0.05,
                                 nodeNames = modNodes,
                                 nodeAttr = moduleColors[inModule])
  
}





