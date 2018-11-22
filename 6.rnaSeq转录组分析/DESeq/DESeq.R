######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")

setwd("C:/Users/y/Desktop/DESeq")
library("DESeq")
data = read.table("input.txt", row.names=1 ,header=T,check.names=F)
design = factor( c( rep("con",2), rep("treat",2)) )
newTab = newCountDataSet( data, design )
newTab = estimateSizeFactors(newTab)

#pdf(file="rawBox.pdf")
boxplot(data,col = "blue",xaxt = "n",outline = F)
newData=counts(newTab, normalized=TRUE )
#pdf(file="normalBox.pdf")
boxplot(newData,col = "red",xaxt = "n",outline = F)

#Working without any replicates
#newTab = estimateDispersions( newTab, method="blind", sharingMode="fit-only",fitType = "local")
#have replicates
newTab = estimateDispersions( newTab, fitType = "local")
diff = nbinomTest( newTab, "con", "treat")
#pdf(file="qValueHist.pdf")
hist(diff$padj, breaks=100, col="skyblue", border="slateblue", main="")
write.table( diff, file="DESeqOut.xls",sep="\t",quote=F,row.names=F)
diffSig = diff[is.na(diff$padj)==FALSE,]
diffSig = diffSig[(diffSig$padj < 0.05 & (diffSig$log2FoldChange>1 | diffSig$log2FoldChange<(-1))),]
diffSig = diffSig[order(diffSig$padj),]
write.table( diffSig, file="DESeqSig.xls",sep="\t",quote=F,row.names=F)

#heatmap
hmExp=log10(newData[diffSig$id,]+0.00001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf")
par(oma=c(3,3,3,5))
heatmap.2(hmMat,col='greenred',trace="none",cexCol=1)

#volcano
pdf(file="vol.pdf")
allDiff=diff[ is.na(diff$padj)==FALSE,]
xMax=max(-log10(allDiff$padj))+1
yMax=10
plot(-log10(allDiff$padj), allDiff$log2FoldChange, xlab="-log10(padj)",ylab="log2FoldChange",main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=subset(allDiff, allDiff$padj<0.05 & allDiff$log2FoldChange>1)
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=subset(allDiff, allDiff$padj<0.05 & allDiff$log2FoldChange<(-1))
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
