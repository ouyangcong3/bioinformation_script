######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\TCGA11_coexpressionMerge")   #设置工作目录（需修改）
miRNA = read.table("diffmiRNAExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)
RNA = read.table("diffmRNAExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)
colnames(miRNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(miRNA))
colnames(RNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(RNA))
sameSample=intersect(colnames(miRNA),colnames(RNA))
merge=rbind(id=sameSample,miRNA[,sameSample],RNA[,sameSample])
write.table(merge,file="merge.txt",sep="\t",quote=F,col.names=F)

miRNALabel=cbind(rownames(miRNA),"miRNA")
RNALabel=cbind(rownames(RNA),"gene")
nodeLabel=rbind(c("ID","Classify"),miRNALabel,RNALabel)
write.table(nodeLabel,file="nodeLabel.txt",sep="\t",quote=F,col.names=F,row.names=F)

