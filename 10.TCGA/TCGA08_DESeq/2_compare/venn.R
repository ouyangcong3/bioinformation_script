######Video source: https://shop119322454.taobao.com
#install.packages("VennDiagram")

file1="edgeR.txt"
file2="DESeq.txt"

setwd("C:\\Users\\YunGu\\Desktop\\compare")
list1=list()
library(VennDiagram)
rt1=read.table(file1,header=F)
list1[["edgeR"]]=as.vector(rt1[,1])
rt2=read.table(file2,header=F)
list1[["DESeq"]]=as.vector(rt2[,1])
venn.diagram(list1,filename="venn.png",fill=c("green","red"),cat.cex=0.6)

sect1=intersect(list1[[1]],list1[[2]])
write.table(file="intersect.xls",sect1,sep="\t",quote=F,col.names=F,row.names=F)
