######Video source: https://shop119322454.taobao.com
dat <- read.table("goAnnResult.txt", header=T,sep="\t",check.names=F,quote = "")
pValue<-numeric(0)
for (i in 1:length(dat[,1])){
    d <- data.frame(allGene=c(dat[i,]$X,dat[i,]$N), diffGene=c(dat[i,]$x,dat[i,]$n))
    a=fisher.test(d)
    pValue=c(pValue,a$p.value)}
qValue<-p.adjust(pValue,method="fdr")
dat2=cbind(dat,pValue,qValue)
dat2=dat2[order(dat2$pValue),]
write.table(dat2,file="goEnrichmentResult.xls",sep="\t",quote=F,row.names = F)

