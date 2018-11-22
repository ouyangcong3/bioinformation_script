######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\GDSC9_corPlot")
input="input.txt"
pThreshold=0.01

rt=read.table(input,sep="\t",header=T,check.names=F)
yMax=ifelse(min(rt$pValue)==0,50,max(-log10(rt$pValue))+1)
tiff(file="cor_pvalue.tiff",width = 14,height = 12,units ="cm",compression="lzw",bg="white",res=1000)
plot(rt$Cor, -log10(rt$pValue),xlab="Pearson correlation coefficient",ylab="-log10(pValue)",
  pch=20, cex=0.4, xlim=c(-1,1),ylim=c(0,yMax))
sig=rt[rt$pValue<pThreshold & rt$Cor>0,]
points(sig$Cor, -log10(sig$pValue), pch=20, col="red",cex=0.4)
sig=rt[rt$pValue<pThreshold & rt$Cor<0,]
points(sig$Cor, -log10(sig$pValue), pch=20, col="green",cex=0.4)
abline(v=0,lty=2,lwd=3)
text(-0.6,2.6,"AT-7519",col="green")
text(-0.7,2.3,"LY317615",col="green")
dev.off()
