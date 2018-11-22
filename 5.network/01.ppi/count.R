######Video source: https://shop119322454.taobao.com
setwd("C:\\Users\\YunGu\\Desktop\\network\\01.ppi")  #dir
rt=read.table("input.txt", header=T, sep="\t", comment.char = "", check.names =FALSE)
tb=table(c(as.vector(rt[,1]),as.vector(rt[,2])))
tb=sort(tb,decreasing =T)
write.table(tb,file="count.xls",sep="\t",quote=F,col.names=F)

n=as.matrix(tb)[1:30,]   #gene number
pdf(file="barplot.pdf")
par(mar=c(3,10,3,3),xpd=T)
bar=barplot(n,horiz=TRUE,col="skyblue",names=FALSE)
text(x=n-1,y=bar,n)
text(x=-0.2,y=bar,label=names(n),xpd=T,pos=2)
dev.off()
