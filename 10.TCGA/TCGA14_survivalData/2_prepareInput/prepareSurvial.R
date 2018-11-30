######Video source: https://shop119322454.taobao.com
#install.packages("hash")

gene='BRCA1'   #基因名（需修改）
clinicalFile="time.txt"
expFile="normalizeExp.txt"

library(hash)
setwd("C:\\Users\\YunGu\\Desktop\\TCGA14_survivalData\\2_prepareInput")     #工作目录（需修改）
rt=read.table(clinicalFile,header=T,check.names=F,sep="\t")
h = hash(keys = rt$id, values = paste(rt$futime,rt$fustat,sep="\t"))

exp=read.table(expFile,header=T,check.names=F,row.names=1,sep="\t")
geneExp=t(exp[gene,])
write.table("sample\tfutime\tfustat\texpression",file="survivalInput.txt",sep="\t"
  ,quote=F,row.names=F,col.names=F)
for(i in rownames(geneExp)){
  j=unlist(strsplit(i,"\\-"))
  if(grepl("^0",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,h[[name3]],geneExp[i,],sep="\t"),file="survivalInput.txt",sep="\t",
                  quote=F,append=TRUE,row.names=F,col.names=F)
                  }
    }
}
