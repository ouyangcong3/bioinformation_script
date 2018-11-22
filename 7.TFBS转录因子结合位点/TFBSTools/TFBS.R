######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#source("https://bioconductor.org/biocLite.R")
#biocLite("TFBSTools")
#biocLite("JASPAR2016")

setwd("C:\\Users\\YunGu\\Desktop\\TFBSTools")   #设置工作目录
fastaFile="gene.fa"                             #fastq文件名

library("TFBSTools")
library(JASPAR2016)
library(Biostrings)

opts = list()
opts[["species"]] <- 9606                      #species id
opts[["matrixtype"]] <- "PWM"
PWMatrixList = getMatrixSet(JASPAR2016, opts)

dnaSet=readDNAStringSet(fastaFile, format="fasta")
dnaTab=as.data.frame(dnaSet)
j=0
for(i in rownames(dnaTab))
{
 j=j+1
 seqName=i
 seqFas=as.character(dnaTab[seqName,])
 subject = DNAString(seqFas)
 sitesetList = searchSeq(PWMatrixList, subject, seqname=seqName, min.score="80%", strand="*")
 if(j==1){
   write.table(writeGFF3(sitesetList,scoreType="relative"),file="result.gff3",quote=F,
   sep="\t",row.names=F)}
 else{
   write.table(writeGFF3(sitesetList,scoreType="relative"),file="result.gff3",quote=F,
   sep="\t",append = TRUE,col.names = FALSE,row.names=F)}
}

