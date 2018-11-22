######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("WGCNA")

library(WGCNA)
setwd("C:\\Users\\YunGu\\Desktop\\network\\03.WGCNA")  #dir

######input######
rt=read.table("input.txt",sep="\t",row.names=1,header=T,check.names=F,quote="!")
datSummary=rownames(rt)
datExpr = t(rt)

######module detection######
ADJ= adjacency(datExpr)
dissTOM=TOMdist(ADJ)
hierTOM = hclust(as.dist(dissTOM),method="average")
colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.8, minSize = 10)
pdf(file="module.pdf")
par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub="")
plotColorUnderTree(hierTOM,colors=data.frame(module=colorh1))
title("Module (branch) color")
dev.off()

######output module######
datME=moduleEigengenes(datExpr,colorh1)[[1]]
color1=rep("grey",dim(datExpr)[[2]])
color1=as.character(colorh1)
datKME=signedKME(datExpr, datME)
datout=data.frame(datSummary, colorNEW=color1,datKME )
write.table(datout, file="OutputModule.xls", sep="\t", row.names=F,quote=F)

######network visual######
exportNetworkToCytoscape(ADJ,edgeFile="edge.txt",nodeFile="node.txt",threshold = 0.9)


########################################################################################


######Relating modules to physiological traits######
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,colorh1)$eigengenes
MEsFemale = orderMEs(MEs0)
datTraits = read.table("survial.txt",sep="\t",row.names="id",header=T,check.names=F,quote="!")
modTraitCor = cor(MEsFemale, datTraits, use = "p")
write.table(file="modPhysiological.cor.xls",modTraitCor,sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(file="modPhysiological.p.xls",modTraitP,sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), " (", signif(modTraitP, 1), ")",sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf(file="modPhysiologicalTraits.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), 
  yLabels = names(MEsFemale),ySymbols = names(MEsFemale), colorLabels = FALSE, 
  colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, 
  cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

######Relating genes to physiological traits######
brownExp=datExpr
brownTraitCor = cor(brownExp, datTraits, use = "p")
write.table(file="genePhysiological.cor.xls",brownTraitCor,sep="\t",quote=F)
brownTraitP = corPvalueStudent(brownTraitCor, nSamples)
write.table(file="genePhysiological.p.xls",brownTraitP,sep="\t",quote=F)
