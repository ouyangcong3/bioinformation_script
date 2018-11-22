#install.packages("RSQLite")
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

setwd("C:\\Users\\guangjun\\Desktop\\symbol2id")

library("org.Hs.eg.db")
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

######ÊÓÆµ×ÊÁÏÍøÖ·£ºwww.xixibio.com
######ÌÔ±¦µêÆÌÍøÖ·£ºhttps://shop119322454.taobao.com
