outputDir <- function(path){
  
  ########################################
  dirRes <- paste0(path,"/Results/")
  if(!file.exists(dirRes)){
    dir.create(dirRes)
  }
  ########################################
  # output 
  
  dirTxt <- paste0(dirRes,"txtFile/")
  if(!file.exists(dirTxt)){
    dir.create(dirTxt)
  }
  
  dirFigure <- paste0(dirRes,"figure/")
  if(!file.exists(dirFigure)){
    dir.create(dirFigure)
  }
  
  dirRdata <- paste0(dirRes,"Rdata/")
  if(!file.exists(dirRdata)){
    dir.create(dirRdata)
  }
  ########################################

  output_directory <- list(dirTxt = dirTxt,
                           dirFigure = dirFigure,
                           dirRdata = dirRdata)
  
  return(output_directory)
}
