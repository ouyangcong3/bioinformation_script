######Video source: https://shop119322454.taobao.com
#install.packages("pheatmap")

setwd("C:\\Users\\YunGu\\Desktop\\heatmap")
rt=read.table("diffExp.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[1:nrow(rt),]   #选择所有的差异基因画热图
#rt=rt[1:50,]        #选择前50个基因
library(pheatmap)
annotation=read.table("group.txt",sep="\t",header=T,row.names=1)

tiff(file="heatmap1.tiff",width = 15,height = 35,units ="cm",compression="lzw",bg="white",res=500)
pheatmap(rt, annotation=annotation,fontsize_row=7,fontsize_col=10)
dev.off()

tiff(file="heatmap2.tiff",width = 15,height = 35,units ="cm",compression="lzw",bg="white",res=500)
pheatmap(rt, annotation=annotation, color = colorRampPalette(c("green", "black", "red"))(50),fontsize_row=8,fontsize_col=10)
dev.off()

tiff(file="heatmap3.tiff",width = 15,height = 35,units ="cm",compression="lzw",bg="white",res=500)
pheatmap(rt, annotation=annotation, display_numbers = TRUE,fontsize_row=8,fontsize_col=10)
dev.off()

tiff(file="heatmap4.tiff",width = 15,height = 35,units ="cm",compression="lzw",bg="white",res=500)
pheatmap(rt, annotation=annotation, cluster_cols = FALSE,fontsize_row=8,fontsize_col=10)
dev.off()

#the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D",
# "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
tiff(file="heatmap5.tiff",width = 15,height = 35,units ="cm",compression="lzw",bg="white",res=500)
pheatmap(rt, annotation=annotation, clustering_method = "median",fontsize_row=8,fontsize_col=10)
dev.off()
