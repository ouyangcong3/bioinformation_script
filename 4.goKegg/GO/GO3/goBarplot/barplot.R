setwd("C:\\Users\\guangjun\\Desktop\\goBarplot")
r1=read.table("bp.txt",sep="\t",header=F,as.is=T,quote="!")
r2=read.table("cc.txt",sep="\t",header=F,as.is=T,quote="!")
r3=read.table("mf.txt",sep="\t",header=F,as.is=T,quote="!")
x1=nrow(r1); x2=nrow(r2); x3=nrow(r3); x=x1+x2+x3
m=c(r1[,2],r2[,2],r3[,2])
l=c(r1[,1],r2[,1],r3[,1])
y=ceiling(max(m)/10)*15
tiff(file="goarplot.tiff",width = 30,height = 20,units ="cm",compression="lzw",bg="white",res=300)
par(mar=c(8,4,3,3),mgp=c(0.8,0.3,0),cex.axis=.7)
barplot(m,beside=T,col=c(rep(rgb(153/255,216/255,201/255),x1),rep(rgb(44/255,127/255,184/255),x2),rep(rgb(201/255,148/255,199/255),x3)),space=0,xaxs='i',yaxs='i',yaxt='n',las=2,names.arg=l,ylab="target genes")
abline(h=0)
par(xpd=T)
lx=max(nchar(l))
y1=5/100*y
y2=10/100*y
segments(0,-y1,0,-y2); segments(0,-y2,x,-y2); segments(x1,-y1,x1,-y2); segments(x1+x2,-y1,x1+x2,-y2); segments(x,-y1,x,-y2)
text(x1/2,-(y2-2.5/10*y2),labels='biological_process',pos=1,cex=1,col=rgb(153/255,216/255,201/255))
text(x1+x2/2,-(y2-2.5/10*y2),labels='cellular_component',pos=1,cex=1,col=rgb(44/255,127/255,184/255))
text(x1+x2+x3/2,-(y2-2.5/10*y2),labels='molecular_function',pos=1,cex=1,col=rgb(201/255,148/255,199/255))
axis(2)
dev.off()

###### ”∆µ◊ ¡œÕ¯÷∑£∫www.xixibio.com
######Ã‘±¶µÍ∆ÃÕ¯÷∑£∫https://shop119322454.taobao.com
