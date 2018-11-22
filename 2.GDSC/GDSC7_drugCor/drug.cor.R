######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\GDSC7_drugCor")
drugFile="drug.txt"
expFile="expression.txt"
drug="1498"

expData=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
drugData=read.table(drugFile,sep="\t",header=T,check.names=F,row.names=1)

outputFile=paste(drug,".cor.xls",sep="")
outTab=data.frame()

for(gene in rownames(expData)){
   corTab=rbind(expData[gene,],drugData[drug,])
   corTab=corTab[,apply(corTab,2,min)>-100]
   x=as.numeric(corTab[gene,])
   y=as.numeric(corTab[drug,])
   Cor=cor(x,y)
   Cor=round(Cor,3)
   outTab=rbind(outTab,cbind(drug,gene,Cor))
  }

corVal=as.numeric(as.vector(outTab[,3]))
corMean=mean(corVal)
corSd=sd(corVal)
zVal=(corVal-corMean)/corSd
pValue=2*pnorm(-abs(zVal))
wt=cbind(outTab,corMean,corSd,zVal,pValue)
wt=wt[order(wt$pValue),]
write.table(file=outputFile,wt,sep="\t",quote=F,row.names=F)
