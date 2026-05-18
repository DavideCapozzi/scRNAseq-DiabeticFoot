relateModstoExt <- function(){
  
  ########################################
  # input parameters
  
  trait_name <- input_parameter$trait_name
  
  module <- input_parameter$module

  filename_dataInput <- output_file$filename_dataInput
  
  filename_networkConstruction <- output_file$filename_networkConstruction
  
  filename_driverGenes <- output_file$filename_driverGenes
  
  filename_geneInfo <- output_file$filename_geneInfo
  ########################################
  
  lnames = load(file = filename_dataInput)
  
  lnames = load(file = filename_networkConstruction)
  ########################################
  
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  ########################################
  #base_path <- "G:/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V1/Results_V1/txtFile/" 
  base_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/sc-FEDE_DAVIDE/FootPlasma_data/analysis_joined/proteomics_analysis/code/project/CeciliaMorgantini/dataset/DiabeticFoot_V3/Results_V3/txtFile/"
  
  # 2. Salva la matrice delle Correlazioni
  write.table(moduleTraitCor, 
              file = paste0(base_path, "moduleTraitCorV3.txt"), 
              sep = "\t", quote = FALSE, col.names = NA)
  
  # 3. Salva la matrice dei P-values
  write.table(moduleTraitPvalue, 
              file = paste0(base_path, "moduleTraitPvalueV3.txt"), 
              sep = "\t", quote = FALSE, col.names = NA)
  
  
  ########################################
  
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)

  par(mar = c(6, 8.5, 3, 3));
  p <- labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))

  ########################################
  
  trait = as.data.frame(datTraits[, trait_name])
  names(trait) = trait_name

  modNames = substring(names(MEs), 3)
  
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(datExpr, trait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  names(geneTraitSignificance) = paste("GS.", names(trait), sep="");
  names(GSPvalue) = paste("p.GS.", names(trait), sep="");
  ########################################
  
  module = module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  pdf(filename_driverGenes,width=6, height=5)
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", trait_name),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  ########################################

  geneInfo0 = data.frame(moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  
  modOrder = order(-abs(cor(MEs, trait, use = "p")));

  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[,2]));
  geneInfo = geneInfo0[geneOrder, ]
  ########################################
  
  write.table(geneInfo, file = filename_geneInfo, sep = "\t", row.names = T, col.names = NA, quote = F)
  
  
}





