######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("WGCNA")

setwd("C:\\Users\\YunGu\\Desktop\\TCGA12_WGCNA\\WGCNA")     #设置工作目录（需修改）
library(WGCNA)
rt=read.table("merge.txt",sep="\t",row.names=1,header=T,check.names=F,quote="!")
datExpr = t(rt)

######select beta value######
powers1=c(seq(1,10,by=1),seq(12,30,by=2))
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers1)[[2]]
cex1=0.7
#pdf(file="softThresholding.pdf")
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")

beta1=9
ADJ= adjacency(datExpr,power=beta1)
vis=exportNetworkToCytoscape(ADJ,edgeFile="edge.txt",nodeFile="node.txt",threshold = 0.8)
