#Vikas Bansal April 2020
#Calculate the double positive cells for a combination of genes (at least one count for each).

set.seed(786)
setwd("/data/vikas/COV_ace/")


.libPaths( c( "/home/vikas/Rlib/", .libPaths()) )
library(Seurat)
library(gplots)
library(reshape2)



pathto.outPlots <- "/data/vikas/COV_ace/OutputPlots2/"
pathto.outData <- "/data/vikas/COV_ace/OutputData/"
outName <- "TestisNumberGeneComb"


Adult_Testis <- readRDS(file = "OutputData/Online_put2/Adult_Testis.rds")

iPSCsDopa.integratedv2 <- Adult_Testis
DefaultAssay(iPSCsDopa.integratedv2) <- "RNA"



GeneReceptors <- c("ACE2","DPP4","BSG","ANPEP","CD209","CLEC4G","CLEC4M")
GeneProtease <- c("TMPRSS2","CTSB","CTSL","FURIN","TMPRSS4","TMPRSS11a","TMPRSS11b")



GeneName <- expand.grid(GeneReceptors = GeneReceptors, GeneProtease = GeneProtease)
GeneName$GeneReceptors <- as.character(GeneName$GeneReceptors)
GeneName$GeneProtease <- as.character(GeneName$GeneProtease)

#ExpressionThresh <- 0


for (i in 1:nrow(GeneName)){
  possibleError <- tryCatch(
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[c(GeneName[i,1],GeneName[i,2]),],
    error=function(e) e
  )
  
  if(!inherits(possibleError, "error")){

    
    
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i,1],]
    Take_cells1[Take_cells1 > 0]=1
    
    Take_cells2 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i,2],]
    Take_cells2[Take_cells2 > 0]=1
    Take_cells <- Take_cells1 + Take_cells2 
    
    Take_PoscellNames <- (names(which(Take_cells>1)))
  }else(Take_PoscellNames <- list())
  
  
  
  Total_cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[,"NewIDsComb"])))
  colnames(Total_cells)[2] <- "TotalNumberofCells"
  
  
  if(length(Take_PoscellNames)>0){Positive_Cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[Take_PoscellNames,"NewIDsComb"])))
  colnames(Positive_Cells)[2] <- "PositiveNumberofCells"}else{
    Positive_Cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[,"NewIDsComb"])))
    colnames(Positive_Cells)[2] <- "PositiveNumberofCells"
    Positive_Cells[,2] <- 0
  }
  
  Merged_positiveCells <- merge(Positive_Cells,Total_cells, by="Var1",all=T)
  Merged_positiveCells[is.na(Merged_positiveCells)] = 0
  
  Merged_positiveCells$PercentageExpressed <- (Merged_positiveCells$PositiveNumberofCells/Merged_positiveCells$TotalNumberofCells)*100
  colnames(Merged_positiveCells)[1] <- "NewIDsComb"
  colnames(Merged_positiveCells)[c(2,4)] <- paste0(colnames(Merged_positiveCells)[c(2,4)],"_",GeneName[i,1],"_",GeneName[i,2])
  
  
  if(i==1){Merged_positiveCellsV2 <- Merged_positiveCells}else{
    Merged_positiveCellsV2 <- merge(Merged_positiveCellsV2,Merged_positiveCells)
  } 
  
}



write.table(Merged_positiveCellsV2,file = paste0(pathto.outData,outName,"_PositiveCells.txt"), row.names=F, quote=F, sep="\t")
