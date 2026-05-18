outputFiles <- function(){
  
  ########################################
  # input parameters
  
  path <- input_parameter$path

  module <- input_parameter$module
  
  output_dir <- outputDir(path)
  ########################################
  # output directories
  
  dirTxt <- output_dir$dirTxt
  dirFigure <- output_dir$dirFigure
  dirRdata <- output_dir$dirRdata
  ########################################
  # output summary file
  
  filename_parameters <- paste0(path,"/parameter.RData")
  ########################################
  # output txtFile
  
  filename_geneInfo <- paste0(dirTxt,"geneInfo.txt")
  filename_cytoscapeInput_edges <- paste0(dirTxt,"cytoscapeInput_edges.txt")
  filename_cytoscapeInput_nodes <- paste0(dirTxt,"cytoscapeInput_nodes.txt")
  ########################################
  # output figure
  
  filename_moduleSize_barplot <- paste0(dirFigure,"moduleSize_barplot.pdf")
  filename_driverGenes <- paste0(dirFigure,"driverGenes_", module, ".pdf")
  ########################################
  # output Rdata
  
  filename_dataInput <- paste0(dirRdata,"dataInput.RData")
  filename_networkConstruction <- paste0(dirRdata,"networkConstruction.RData")
  ########################################

  output_file <- list(filename_geneInfo = filename_geneInfo,
                      filename_cytoscapeInput_edges = filename_cytoscapeInput_edges,
                      filename_cytoscapeInput_nodes = filename_cytoscapeInput_nodes,
                      filename_moduleSize_barplot = filename_moduleSize_barplot,
                      filename_driverGenes = filename_driverGenes,
                      filename_dataInput = filename_dataInput,
                      filename_networkConstruction = filename_networkConstruction)
  
  return(output_file)
}



