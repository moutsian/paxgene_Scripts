#PCA for RNA-seq Vaccines pilot and GTEx
#Author: Javier Gutierrez-Achury - WTSI 
#usage: R --vanilla --args matrixInGTEx "ExpressionMatrixGTEx" matrixInPilot "ExpressionMatrixPilot" outputPcs "ExpressionMatrixPCS" < /lustre/scratch113/teams/anderson/users/jga/ToolBox/ScriptsR/PCA.ExpressionData.R > PCA.ExpressionData.Rout

### MDS in GTEx + PilotData

library('batch')
library(data.table)
library(MASS)
library(ggplot2)

matrixInGTEx = " "
matrixInPilot = " "
outputPcs = " "

parseCommandArgs()
rpkm2tpm = function(data){
  if(class(data[,1]) == "character"){
    geneNames = data[,1]
    colnames = data
    dataToWork = data[,c(-1)]
  } else {
    dataToWork = data
  }
  numberRows = dim(dataToWork)[1]
  numberCols = dim(dataToWork)[2]
  matrixToFill = matrix(c(0), nrow=numberRows , ncol=numberCols)
  vectorOfSums = colSums(dataToWork)
  for(cols in 1:numberCols){
    matrixToFill[,cols] = dataToWork[,cols]/vectorOfSums[cols]*10^6
  }

  if(class(data[,1]) == "character"){
    IDs = 
    rownames(matrixToFill) = data[,1]
    colnames(matrixToFill) = colnames(data)[c(-1)]
  }

  return(matrixToFill)
}


#This is RPKM data
gtexData.rpkm = fread(matrixInGTEx, header=TRUE)
gtexData.rpkm.df = as.data.frame(gtexData.rpkm)
gtexData.rpkm.df = gtexData.rpkm.df[,c(-2)] # remove Ensembl genes IDs gtexData.rpkm.df[,c(-1)]- If want to keep Ensembl, gtexData.rpkm.df[,c(-2)]
rownames(gtexData.rpkm.df) = gtexData.rpkm.df$Name
gtexData.rpkm.df = gtexData.rpkm.df[,c(-1)]
gtexData.rpkm.df.filtered = gtexData.rpkm.df[rowSums(gtexData.rpkm.df) >= 0.05,]
print(dim(gtexData.rpkm.df.filtered))
print("GTEx data in")
#print(head(gtexData.rpkm.df.filtered))
# CHECKING THIS PART!!!! if using RPKM directly for the rest of the analysis, uncomment the next two lines, otherwise... comment!
#rownames(gtexData.rpkm.df) = gtexData.rpkm.df$Description
#gtexData.rpkm.df = gtexData.rpkm.df[,c(-1)] # remove Ensembl genes IDs

## convert RPKM to TPM if Pilot-Vaccines data is TPM. If pilot-Vaccines is already RPKM = DO NOT USE
#gtexData.tpm = rpkm2tpm(gtexData.rpkm.df)


# This is already RPKM
pilotData = fread(matrixInPilot, header=T)
pilotData.df = as.data.frame(pilotData)
rownames(pilotData.df) = pilotData.df$Name
pilotData.df = pilotData.df[,c(-1)]
print(dim(pilotData.df))
print("Pilot data in")
#print(head(pilotData.df))

#rownames(pilotData.df) = pilotData.df$Names
#pilotData.df = pilotData.df[,c(-1)]

#if using RPKM
gtexData.PilotData.merged = merge(gtexData.rpkm.df.filtered, pilotData.df, by="row.names")
print("Pilot data merged")
#if using TPM
#gtexData.PilotData.merged.filtered = gtexData.PilotData.merged[rowSums(gtexData.PilotData.merged) > 0.05,] # this filter - 0.05 - can change

#write.table(gtexData.PilotData.merged, file=outputPcs, col.names=T, row.names=F, quote=F, sep="\t")

rownames(gtexData.PilotData.merged) = gtexData.PilotData.merged$Row.names
gtexData.PilotData.merged  = gtexData.PilotData.merged[,c(-1)]
gtexData.PilotData.merged.log2 = log2(gtexData.PilotData.merged+1)
gtexData.PilotData.merged.log2.t = t(gtexData.PilotData.merged.log2)

gtexData.PilotData.merged.log2.dist = dist(gtexData.PilotData.merged.log2.t)
gtexData.PilotData.merged.log2.dist.msd = isoMDS(gtexData.PilotData.merged.log2.dist, k=10)

write.table(gtexData.PilotData.merged.log2.dist.msd$points, file=outputPcs, col.names=T, row.names=T, quote=F, sep="\t")



