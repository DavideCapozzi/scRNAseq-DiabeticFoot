inputFiles <- function(){
  
  ########################################
  # input parameters
  
  filename_data <- input_parameter$filename_data
  filename_traits <- input_parameter$filename_traits
  ########################################
  # matrix
  
  data <- read.table(filename_data, header = T, row.names = 1, sep = '\t', check.names = F, quote = "")
  
  traits <- read.table(filename_traits, header = T, row.names = 1, sep = '\t', check.names = F, quote = "")

  ########################################
  
  input_file <- list(data = data,
                     traits = traits)
  
  return(input_file)
}
