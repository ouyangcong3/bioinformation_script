######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\network\\02.cor")  #dir
rt=read.table("data.txt",sep="\t",header=T,check.names=F,row.names=1)
i=1  #gene 1
j=2  #gene 2 
x=as.numeric(rt[i,])
y=as.numeric(rt[j,])
xName=row.names(rt[i,])
yName=row.names(rt[j,])
R=cor(x,y)
R=round(R,3)

pdf(file="cor.pdf")
z=lm(y~x)
plot(x,y, type="p",pch=16,cex=0.8, col="blue",main=paste("Pearson's correlation=",R,sep=""),
xlab=paste(xName,"mRNA expression"),ylab=paste(yName,"mRNA expression"))
lines(x,fitted(z),col="red")
dev.off()
