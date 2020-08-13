#Vikas Bansal July 2020
#Calculate the significance (Fisher's exact test) number/percentage of positive cells for a gene (at least one count) and create bonferroni corrected heatmap with star.

set.seed(786)
setwd("/data/vikas/COV_ace/")


.libPaths( c( "/home/vikas/Rlib/", .libPaths()) )
library(Seurat)
library(gplots)
library(reshape2)


pathto.outPlots <- "/data/vikas/COV_ace/OutputPlotsRevision/"
pathto.outData <- "/data/vikas/COV_ace/OutputPlotsRevision/"
outName <- "TestisStar"


Adult_Testis <- readRDS(file = "OutputData/Online_put2/Adult_Testis.rds")

GeneList <- c("TMPRSS2","CD209","CLEC4G","CLEC4M","ACE2","DPP4","BSG","ANPEP","LY6E","IFITM3","IFITM2","IFITM1")


#Can read any seurat object used in the manuscript
iPSCsDopa.integratedv2 <- Adult_Testis 
GeneName <- GeneList



for (i in 1:length(GeneName)){
  possibleError <- tryCatch(
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[c(GeneName[i],GeneName[i]),],
    error=function(e) e
  )
  
  #error handling
  if(!inherits(possibleError, "error")){
   
    
    
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i],]
    Take_cells1[Take_cells1 > 0]=1
    
    Take_cells2 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i],]
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
  colnames(Merged_positiveCells)[c(2,4)] <- paste0(colnames(Merged_positiveCells)[c(2,4)],"_",GeneName[i])
  
  
  if(i==1){Merged_positiveCellsV2 <- Merged_positiveCells}else{
    Merged_positiveCellsV2 <- merge(Merged_positiveCellsV2,Merged_positiveCells)
  } 
  
}




Merged_positiveCellsV3 <- Merged_positiveCellsV2[,-grep("Percentage",colnames(Merged_positiveCellsV2))]

colnames(Merged_positiveCellsV3) <- gsub("PositiveNumberofCells_","",colnames(Merged_positiveCellsV3)) 


Merged_positiveCellsHeatmap <- Merged_positiveCellsV3

for(cellTypeNum in 1:nrow(Merged_positiveCellsV3)){
  
  for(i in 1:length(GeneName)){
    fish_stat <- fisher.test(matrix(c(Merged_positiveCellsV3[cellTypeNum,GeneName[i]],Merged_positiveCellsV3[cellTypeNum,"TotalNumberofCells"]-Merged_positiveCellsV3[cellTypeNum,GeneName[i]],
                                      sum(Merged_positiveCellsV3[-cellTypeNum,GeneName[i]]),sum(Merged_positiveCellsV3[-cellTypeNum,"TotalNumberofCells"])-sum(Merged_positiveCellsV3[-cellTypeNum,GeneName[i]])),
                                    nrow=2), alternative="greater")
    Merged_positiveCellsHeatmap[cellTypeNum,GeneName[i]] <- fish_stat$p.value
    
  }
  
}



SigGenesToUseMeltedReshaped <- as.matrix(Merged_positiveCellsHeatmap[,-c(1:2)])
rownames(SigGenesToUseMeltedReshaped) <- Merged_positiveCellsHeatmap[,1]

Pval_mat2=SigGenesToUseMeltedReshaped

#Put star that survived bonferroni correct pvalue
Pval_mat2[(Pval_mat2<=(0.05/(nrow(Pval_mat2)*ncol(Pval_mat2))))]="*"
Pval_mat2[(!Pval_mat2=="*")]=""



pal=colorRampPalette(c("blue","white","grey"))
col_my=pal(18)


pdf(paste0(pathto.outPlots,"Enrichment_Heatmap_SigPercentage_",(0.05/(nrow(Pval_mat2)*ncol(Pval_mat2))),"_",outName,".pdf"), width=10, height=10)


heatmap.2(SigGenesToUseMeltedReshaped,scale="none", cellnote=Pval_mat2, key.xlab = "Fisher's Pvalue", notecex=3, notecol="black", col=col_my,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexCol=1.5,cexRow=1.5, dendrogram="none", Rowv=F, Colv=F, keysize=1,margin=c(12,18))

dev.off()




