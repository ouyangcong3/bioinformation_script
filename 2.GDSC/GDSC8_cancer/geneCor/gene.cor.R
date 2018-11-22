######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\GDSC8_cancer\\geneCor")
drugFile="drug.txt"
expFile="expression.txt"
gene="BRAF"
sampleFile="sample.txt"

expData=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
drugData=read.table(drugFile,sep="\t",header=T,check.names=F,row.names=1)
sample=read.table(sampleFile,header=F,check.names=F)
samSample=intersect(sample[,1],colnames(expData))
expData=expData[,samSample]
drugData=drugData[,samSample]

outputFile=paste(gene,".cor.xls",sep="")
outTab=data.frame()

for(drug in rownames(drugData)){
   corTab=rbind(expData[gene,],drugData[drug,])
   corTab=corTab[,apply(corTab,2,min)>-100]
   x=as.numeric(corTab[gene,])
   y=as.numeric(corTab[drug,])
   Cor=cor(x,y)
   Cor=round(Cor,3)
   outTab=rbind(outTab,cbind(gene,drug,Cor))
  }

corVal=as.numeric(as.vector(outTab[,3]))
corMean=mean(corVal)
corSd=sd(corVal)
zVal=(corVal-corMean)/corSd
pValue=2*pnorm(-abs(zVal))
wt=cbind(outTab,corMean,corSd,zVal,pValue)
wt=wt[order(wt$pValue),]
write.table(file=outputFile,wt,sep="\t",quote=F,row.names=F)
