######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")
#biocLite("Sushi")

library("Sushi")

chrom = "chr21"
chromstart = 45202761
chromend = 45203527
treatFile = "/home/yun/bio/chip/ERR990/sample_MACS_bedGraph/treat/sample_treat_afterfiting_chr21.bdg"
controlFile = "/home/yun/bio/chip/ERR990/sample_MACS_bedGraph/control/sample_control_afterfiting_chr21.bdg"

pdf(file="peakView.pdf",width=12)
par(mar=c(4,5,3,3))

date=read.table(treatFile,sep="\t",header=T,skip=1)
date1=read.table(controlFile,sep="\t",header=T,skip=1)
yMax=max(c(date[,4],date1[,4]))
plotBedgraph(date,chrom,chromstart,chromend,range=c(0,yMax),
			              transparency=.50,color=SushiColors(2)(2)[1])
plotBedgraph(date1,chrom,chromstart,chromend,range=c(0,yMax),
			             transparency=.50,color=SushiColors(2)(2)[2],overlay=TRUE,)

labelgenome(chrom,chromstart,chromend,n=5,scale="bp")
legend("topright",inset=0.025,legend=c("Treat","Control"),
	          fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,
			         cex=1.0)
axis(side=2,las=2,tcl=.2)
mtext("Read Depth",side=2,line=3,cex=1,font=2)
dev.off()

