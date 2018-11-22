#install.packages("colorspace")
#install.packages("stringi")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DOSE")
#biocLite("clusterProfiler")
#biocLite("pathview")

setwd("C:\\Users\\guangjun\\Desktop\\clusterProfiler")
library("clusterProfiler")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

geneFC=rt$logFC
gene=rt$entrezID
names(geneFC)=gene

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="barplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(kk, drop = TRUE, showCategory = 12)
dev.off()

#点图
tiff(file="dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(kk)
dev.off()

#通路图
library("pathview")
keggxls=read.table("KEGG.txt",sep="\t",header=T)
for(i in keggxls$ID){
  pv.out <- pathview(gene.data = geneFC, pathway.id = i, species = "hsa", out.suffix = "pathview")
}

######视频资料网址：www.xixibio.com
######淘宝店铺网址：https://shop119322454.taobao.com
