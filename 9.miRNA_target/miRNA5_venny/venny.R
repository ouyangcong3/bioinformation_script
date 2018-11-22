#video source: https://shop119322454.taobao.com
#install.packages("VennDiagram")

setwd("C:\\Users\\guangjun\\Desktop\\venny")
files=dir()
files=grep("txt",files,value=T)
targetList=list()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=F)
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    targetList[[header[1]]]=as.vector(rt[,1])
    uniqLength=length(unique(as.vector(rt[,1])))
    print(paste(header[1],uniqLength,sep=" "))
}

library(VennDiagram)
venn.diagram(targetList,filename="venny.tiff",imagetype = "tiff",
             fill=rainbow(length(targetList)),cat.cex=1)
intersectGenes=Reduce(intersect,targetList)
write.table(file="target.xls",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
