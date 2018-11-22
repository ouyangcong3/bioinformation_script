#video source: https://shop119322454.taobao.com
#install.packages("ggplot2")
setwd("C:\\Users\\guangjun\\Desktop\\miRNA6_GO")
library(ggplot2)
dat=read.table("input.txt", header = T,sep="\t",check.names=F)
tiff(file="GO.tiff",width = 35,height = 18,units ="cm",compression="lzw",bg="white",res=600)
ggplot(data=dat)+geom_bar(aes(x=Term, y=Count, fill=PValue), stat='identity') + coord_flip() + scale_fill_gradient(low="blue", high = "red",name = "PValue") + xlab("") + ylab("Gene count") + theme(axis.text.x=element_text(color="black", size=12), axis.text.y=element_text(color="black", size=12)) + scale_x_discrete(expand=c(0,0))
dev.off()
