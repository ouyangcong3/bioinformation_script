######Video source: https://shop119322454.taobao.com
dat=read.table("peakRegionStat.xls",sep="\t")
r=round(100*dat[,2]/sum(dat[,2]),2)
if(sum(r)==100){ 
		r=r
}else{
	x=which(r==max(r))[1]
	r[x]=100-sum(r[-x])
}
ratio=paste(r,"%)",sep="")
label=paste(dat[,1],ratio,sep="(")
pdf(file="peakRegionStat.pdf",width=10)
pie(dat[,2], labels = label, edges = 200, radius = 0.8)

