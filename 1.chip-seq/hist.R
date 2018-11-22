######Video source: https://shop119322454.taobao.com
date=read.table("sample_peaks.xls",sep="\t",header=T,skip=23)
nrow(date)
date=date[date$length<=1000,]
nrow(date)
pdf(file="peakHist.pdf")
hist(date$length, breaks = 50, col = "lightblue", xlim=c(0,1000),xlab="peak length", main="Histogram of peak length")
dev.off()

