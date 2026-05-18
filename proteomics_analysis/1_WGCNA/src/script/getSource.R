getSource <- function(){
  
  # input
  source("src/script/input/config.R")
  source("src/script/input/inputFiles.R")
  
  # output
  source("src/script/output/outputDir.R")
  source("src/script/output/outputFiles.R")
  
  # lib
  source("src/script/lib/dataInput.R")
  source("src/script/lib/networkConstruction.R")
  source("src/script/lib/relateModstoExt.R")
  source("src/script/lib/visualization.R")
  source("src/script/lib/exportNetwork.R")
  
}
