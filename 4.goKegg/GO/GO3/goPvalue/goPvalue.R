######Video source: https://shop119322454.taobao.com
#install.packages("ggplot2")

setwd("C:\\Users\\guangjun\\Desktop\\goPvalue")
library(ggplot2)
dat <- read.table("input.txt", header = T,sep="\t")
tiff(file="goPvalue.tiff",width = 30,height = 15,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat)+geom_bar(aes(x=Term, y=Count, fill=-log10(PValue)), stat='identity')+
  coord_flip() + scale_fill_gradient(low="red", high = "blue")+ 
  xlab("") + ylab("") + theme(axis.text.x=element_text(color="black", size=12),
  axis.text.y=element_text(color="black", size=12)) + 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))
dev.off()

######ÊÓÆµ×ÊÁÏÍøÖ·£ºwww.xixibio.com
######ÌÔ±¦µêÆÌÍøÖ·£ºhttps://shop119322454.taobao.com

