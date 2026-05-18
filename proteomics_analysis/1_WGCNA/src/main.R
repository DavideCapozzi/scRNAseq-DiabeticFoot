rm(list = ls())
#in V3 è stata commentata una parte dentro dataInput.R

options(stringsAsFactors = FALSE)

#setwd("G:/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code")
setwd("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code")
######################################
library(WGCNA)
library(ggplot2)

source("src/script/getSource.R")
######################################
getSource()
input_parameter <- config()
input_file <- inputFiles()
output_file <- outputFiles()
######################################

dataInput()

networkConstruction()

relateModstoExt()

visualization()

exportNetwork()
