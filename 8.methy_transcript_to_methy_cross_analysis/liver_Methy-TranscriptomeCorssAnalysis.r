suppressMessages(library(limma))
suppressMessages(library(survival))
suppressMessages(library(hash))
suppressMessages(library(car))

#1.对甲基化矩阵进行标准化及整合基因名#########################################################################################
if(file.exists("normalizeMethyProbes.txt")){	#标准化甲基化矩阵文件是否存在,若存在
normalData<-read.table("normalizeMethyProbes.txt",header=T,check.names=F,row.names=1,sep="\t")
}else{	#若不存在
normalNum=45	#正常样品的数目
tumorNum=414	#癌症样品的数目
outTab=data.frame()
grade=c(rep(1,normalNum),rep(2,tumorNum))
Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
rt=read.table("MethylationMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)	#数据框转换成矩阵
rownames(rt)=rt[,1]	#矩阵行名设为与第一列内容(探针名)相同
exp=rt[,2:ncol(rt)]	#删除第一列
dimnames=list(rownames(exp),colnames(exp))	#设定因子格式
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)	#对内容转换为因子格式
data=avereps(data)	#去除探针重复项,取平均
data=data[rowMeans(data)>0,]	#极大值
data=normalizeBetweenArrays(data)	#标准化
normalData=as.data.frame(cbind(id=row.names(data),data))	#将探针编号合并到标准化后的表格的第一列
Probesann<-read.table("probesann.txt")	#读入探针注释列表,来源于TCGAMehtylation.sh脚本中第27步
normalData$ann<-as.vector(t(Probesann))	#将探针注释信息并入标化后甲基化矩阵的最后一列
rownames(normalData)<-as.vector(read.table("Probesnames.txt")[,1])	#由于标准化后无检测信号的探针的探针名被NA替换,不能行名化,所以重新输入来自TCGAMehtylation.sh脚本中第36步的探针名列表 
normalData<-normalData[,-1]	#去掉第一列(探针名)
write.table(normalData,file="normalizeMethyProbes.txt",sep="\t",quote=F)
}

#2.对标准化后的甲基化矩阵中提取位于TSS附近2kb范围内的所有病例甲基化值
TSSProbes<-read.table("TSS2kbProbes.txt")	#读入位于TSS附近2kb范围内的探针信息,来源于TCGAMehtylation.sh脚本中第45步
normalData<-normalData[as.vector(TSSProbes[,1]),]	#提取对应探针的所有病例甲基化值

#3.制作成基因-甲基化值矩阵#########################################################################################
if(!file.exists("TCGA_liver_normalizeMethy.txt")){	#判断基因-甲基化值文件是否"不"存在
genename<-read.table("genename.txt")	#从其他文件提取
write.table(t(as.data.frame(c("geneID",colnames(normalData[,-ncol(normalData)])))),"TCGA_liver_normalizeMethy.txt",sep="\t",quote=F,row.names=F,col.names=F)	#将第一行标题输出
for(i in as.vector(genename[,1])){	#循环提取基因名
colmean<-colMeans(normalData[grepl(i,normalData[,ncol(normalData)]),][,-ncol(normalData)],na.rm=T)	#计算注释为该基因的所有探针的beta均值
write.table(t(as.data.frame(c(i,colmean))),"TCGA_liver_normalizeMethy.txt",sep="\t",quote=F,append=TRUE,row.names=F,col.names=F)
}
exp=read.table("TCGA_liver_normalizeMethy.txt",header=T,check.names=F,row.names=1,sep="\t")	#读入基因-甲基化值矩阵
exp[is.na(exp)]<-0	#将NA,NaN替换成0
write.table(exp,"TCGA_liver_normalizeMethy.txt",sep="\t",quote=F)	#输出为新的矩阵文件
}

#4.分析文件读入#########################################################################################
Methyfile="TCGA_liver_normalizeMethy.txt"	#标准化后的甲基化文件
mRNAFile="normalizeExp.txt"	#标准化后的转录组文件
clinicalFile="time.txt"	#生存时间文件,由survival_time.pl和.json文件运行产生
exp=read.table(Methyfile,header=T,check.names=F,row.names=1,sep="\t")	#读取甲基化矩阵
mRNAexp=read.table(mRNAFile,header=T,check.names=F,sep="\t",row.names=1)	#读取转录组矩阵
time=read.table(clinicalFile,header=T,check.names=F,sep="\t")	#读取时间文件矩阵

#5.甲基化数据和转录组数据交集#########################################################################################
colnames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(exp))	#将TCGA甲基化病人编号格式改成统一的病人编号格式
colnames(mRNAexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(mRNAexp))	#将TCGA转录组病人编号格式改成统一的病人编号格式
samegenenames=intersect(rownames(mRNAexp),rownames(exp))	#提取交集基因名列表
sameSample=intersect(colnames(exp),colnames(mRNAexp))	#提取交集病例编号列表
normalNum<-length(grep("11A",gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)","\\4",sameSample)))	#甲基化病例和转录组病例交集后的正常病例个数

h = hash(keys = time$id, values = paste(time$futime,time$fustat,sep="\t"))
write.table("sample\tfutime\tfustat",file="survivalInput.txt",sep="\t",quote=F,row.names=F,col.names=F)
for(i in sameSample){
  j=unlist(strsplit(i,"\\-"))
  if(grepl("^0",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,h[[name3]],sep="\t"),file="survivalInput.txt",sep="\t",quote=F,append=TRUE,row.names=F,col.names=F)
                  }
    }
}	#循环输出甲基化-转录组交集癌症病例与生存时间病例的交集病例编号
IntersectpatSample<-as.vector(read.table("survivalInput.txt",header=T,sep="\t")[,1])	#读取甲基化-转录组-生存时间病例编号
sameSample<-c(sameSample[1:normalNum],intersect(IntersectpatSample,sameSample[-(1:normalNum)]))	#合并对照组病例编号和生存时间交集的癌症病例编号
tumorNum<-length(intersect(IntersectpatSample,sameSample[-(1:normalNum)]))	#甲基化病例,转录组病例,生存时间病例交集后的癌症病例个数

Intersectexp<-exp[samegenenames,sameSample]	#提取交集甲基化矩阵
IntersectmRNAexp<-mRNAexp[samegenenames,sameSample]	#提取交集转录组矩阵
rownames(exp)=paste(rownames(exp),"methy",sep="|")	#甲基化矩阵中基因名加methy标签
rownames(mRNAexp)=paste(rownames(mRNAexp),"exp",sep="|")	#转录组矩阵中基因名加exp标签
sametagmethygn<-paste(rownames(Intersectexp),"methy",sep="|")	#为交集基因名加methy标签
sametagmRNAgn<-paste(rownames(IntersectmRNAexp),"exp",sep="|")	#为交集基因名加exp标签
mergematrix=rbind(id=sameSample,exp[sametagmethygn,sameSample],mRNAexp[sametagmRNAgn,sameSample])	#提取交集并合并成同一矩阵
write.table(mergematrix,file="merge.txt",sep="\t",quote=F,col.names=F)

seqNum<-1:length(colnames(exp))	#制作和样品量一致的数列(样品名编号)
sameSampleSeq=t(data.frame(seqNum))	#创建数据框
colnames(sameSampleSeq)=colnames(exp)	#将该数据框行名改成甲基化矩阵的列名(样品名)
IntersectNum<-as.numeric(sameSampleSeq[,sameSample])	#提取甲基化-转录组-生存时间交集的样品对应的序号,为第10步层次聚类提取基因对应甲基化探针提供序号

#6.交集甲基化差异分析#########################################################################################
if(file.exists("allGeneMethyDiff.txt")){	#判断甲基化和转录组交集文件是否存在,若存在
outTab=read.table("allGeneMethyDiff.txt",sep="\t",header=T,check.names=F)	#读取甲基化和转录组交集矩阵
}else{	#若不存在

fdrFilter=0.05	#挑选差异基因的fdr阈值

data<-Intersectexp	#将交集甲基化矩阵作为输入数据
outTab=data.frame()
grade=c(rep(1,normalNum),rep(2,tumorNum))
Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))

for(i in row.names(data)){	#基因名循环进行差异分析
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)

  normalGeneMeans=mean(as.numeric(data[i,1:normalNum]))
  tumorGeneMeans=mean(as.numeric(data[i,(normalNum+1):ncol(data)]))
  logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)  
  pvalue=wilcoxTest$p.value
  normalMed=median(as.numeric(data[i,1:normalNum]))
  tumorMed=median(as.numeric(data[i,(normalNum+1):ncol(data)]))
  diffMed=tumorMed-normalMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,
		               cbind(gene=i,
		               normalMean=normalGeneMeans,
		               TumorMean=tumorGeneMeans,
		               logFC=logFC,
		               pValue=pvalue))
	 }
}

pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")	#对p值进行矫正
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="allGeneMethyDiff.txt",sep="\t",row.names=F,quote=F)	#输出所有基因的甲基化差异情况
}
MethyDiffgenenames<-as.character(outTab$gene)

#7.提取甲基化-转录组相关性小于-0.2的基因#########################################################################################

corthreshold=-0.2

if(file.exists("Correlationresult.txt")){	#相关性矩阵文件是否存在,若存在
corresult<-read.table("Correlationresult.txt",header=T,check.names=F,row.names=1,sep="\t")
}else{	#若不存在
write.table("id\tcor\tpvalue",file="Correlationresult.txt",sep="\t",quote=F,row.names=F,col.names=F)
for (gene in MethyDiffgenenames) {
methyGene=paste(gene,"|methy",sep="")
expGene=paste(gene,"|exp",sep="")
i=mergematrix[methyGene,]
j=mergematrix[expGene,]
if(sum(as.numeric(i))==0 | sum(as.numeric(t(j)[,1]))==0){	#若该基因甲基化值或者转录组表达值其中一种缺失
cor<-pval<-NA
}else{
x=as.numeric(i)
y=log2(as.numeric(j)+1)
corT=cor.test(x,y)
methyGeneName=unlist(strsplit(methyGene,"\\|",))[1]
expGeneName=unlist(strsplit(expGene,"\\|",))[1]
cor=corT$estimate
cor=round(cor,3)
pvalue=corT$p.value
pval=signif(pvalue,4)
pval=format(pval, scientific = TRUE)
}
write.table(paste(gene,cor,pval,sep="\t"),file="Correlationresult.txt",row.names=F,col.names=F,append=TRUE,quote=F)
}
corresult<-read.table("Correlationresult.txt",header=T,check.names=F,row.names=1,sep="\t")
corresult<-na.omit(corresult)	#去除NA行
}

filtergenenames<-rownames(corresult[as.numeric(as.vector(corresult$cor))<=corthreshold & as.numeric(as.vector(corresult$pvalue))<0.05,])

for (gene in filtergenenames) {
cor<-as.vector(corresult[gene,1])
pval<-as.vector(corresult[gene,2])

#8.交集转录组差异表达分析#########################################################################################
group.exp.inter<-as.numeric(as.vector(IntersectmRNAexp[gene,]))	#提取该基因的表达量向量
nor.inter=group.exp.inter[1:normalNum]	#提取交集转录组正常组病人表达向量
pat.inter=group.exp.inter[-(1:normalNum)]	#提取交集转录组癌症病人表达向量
factorgroup<-factor(c(rep(1,normalNum),rep(2,tumorNum)))	#创建因子进行方差齐性检验
if(leveneTest(group.exp.inter~factorgroup)[[3]][1]>0.05){	#levene正态性检验,若符合正态性
RNAfcpv<-t.test(nor.inter,pat.inter,alternative = "two.sided")[[3]]	#student's t检验
}else{	#若不符合正态性
rand<-rnorm(length(group.exp.inter),sd=1e-6)	#即避免数据重复导致的wilcox.test"无法精确计算带连结的p值"错误,构建等宽随机数向量
group.exp.inter<-group.exp.inter+rand	#矩阵加上随机数向量
RNAfcpv<-wilcox.test(group.exp.inter[1:normalNum],group.exp.inter[-(1:normalNum)])[[3]]	#wilcox秩和检验的p值
}
RNAlog2fc<-log(mean(pat.inter)/mean(nor.inter),2)	#计算差异倍数

#9.层次聚类分析#########################################################################################
matrix<-normalData[grepl(gene,as.vector(normalData$ann)),][,IntersectNum]	#提取基因对应甲基化探针
a<-as.data.frame(t(matrix))	#转置矩阵,每行为甲基化病例
rownames(a)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",rownames(a))	#将TCGA甲基化病人编号格式改成符合生存分析输入文件中甲基化病人编号的格式
a0<-cbind(id=rownames(a),a)	#将病例号插入矩阵第一列
a1<-a0[-(1:normalNum),]	#去除正常对照组甲基化病人数据
rownames(a1)=NULL	#重新定义行名
out.dist=dist(a1[,-1],method="euclidean")	#对去除正常对照组甲基化病人的数据进行欧氏距离计算 
out.hclust=hclust(out.dist)	#层次聚类分析
out.id=cutree(out.hclust,k=3)	#按定义k组将聚类结果分组
a2<-a1[out.hclust["order"][[1]],]	#按聚类排列
a1$state<-out.id	#将聚类分组编号写入表格
a1<-a1[order(a1$state),]	#表格按聚类分组编号排列
a1$mean<-rowMeans(a1[2:(ncol(a1)-1)])	#将每个探针的平均值输入到a1最后一列,并命名为mean


if(nrow(a1[a1$state==2,])<=1|nrow(a1[a1$state==3,])<=1){
num1<-num2<-num3<-pValue1<-pValue2<-group1log2fc<-group1fcpv<-group2log2fc<-group2fcpv<-group3log2fc<-group3fcpv<-NA
}else{	#在原聚类分组不变的情况下对聚类分组编号重排,将甲基化程度最接近正常组的分为第一组,差异最大分为第三组
g1<-rbind(1,mean(a1[a1$state==1,ncol(a1)]))
g2<-rbind(2,mean(a1[a1$state==2,ncol(a1)]))
g3<-rbind(3,mean(a1[a1$state==3,ncol(a1)]))
sortmat<-data.frame(g1,g2,g3)
if(RNAlog2fc<0){
sortmat<-sortmat[,order(sortmat[2,])]
}else{
sortmat<-sortmat[,order(sortmat[2,],decreasing=T)]
}
sortmat[3,]<-c("A","B","C")
a1[a1$state==sortmat[1,1],ncol(a1)-1]<-rep(sortmat[3,1],nrow(a1[a1$state==sortmat[1,1],]))
a1[a1$state==sortmat[1,2],ncol(a1)-1]<-rep(sortmat[3,2],nrow(a1[a1$state==sortmat[1,2],]))
a1[a1$state==sortmat[1,3],ncol(a1)-1]<-rep(sortmat[3,3],nrow(a1[a1$state==sortmat[1,3],]))

a1[a1$state==sortmat[3,1],ncol(a1)-1]<-rep(1,nrow(a1[a1$state==sortmat[3,1],]))
a1[a1$state==sortmat[3,2],ncol(a1)-1]<-rep(2,nrow(a1[a1$state==sortmat[3,2],]))
a1[a1$state==sortmat[3,3],ncol(a1)-1]<-rep(3,nrow(a1[a1$state==sortmat[3,3],]))

a3<-a1[,c(1,ncol(a1)-1)]	#构建甲基化病人编号和分组编号矩阵
colnames(a3)<-c("sample","state")	#将列名改成sample和state方便下面整合生存时间列表

#10.聚类后各组生存分析#########################################################################################
rt=read.table("survivalInput.txt",header=T,sep="\t")
rt$futime=rt$futime/365	#如果以月为单位，除以30；以年为单位，除以365
rt<-merge(rt,a3,by = intersect(names(rt)[1],names(a3)[1]))	#将聚类分组编号整合至对应病人编号的生存时间列表上
num1<-nrow(rt[rt$state==1,])	#第一类个数
num2<-nrow(rt[rt$state==2,])	#第二类个数
num3<-nrow(rt[rt$state==3,])	#第三类个数
rtg1<-rt[rt$state!=3,]	#选择聚类后的第一组和第二组的生存时间数据
group1<-rtg1$state==2	#提取第二组的生存时间数据(作为生存分析函数的输入数据)
diff1=survdiff(Surv(futime, fustat) ~group1,data = rtg1)	#第一组和第二组的生存分析
pValue1=1-pchisq(diff1$chisq,df=1)	#提取前两组病人生存分析的p值
pValue1=round(pValue1,3)	#设置p值的输出格式
#fit1=survfit(Surv(futime, fustat) ~ group1, data = rtg1)
#summary(fit1)	#查看五年生存率
rtg2<-rt[rt$state!=2,]	#选择聚类后的第一组和第三组的生存时间数据
group2<-rtg2$state==3	#提取第三组的生存时间数据(作为生存分析函数的输入数据)
diff2=survdiff(Surv(futime, fustat) ~group2,data = rtg2)	#第一组和第三组的生存分析
pValue2=1-pchisq(diff2$chisq,df=1)	#提取第一组和第三组生存分析的p值
pValue2=round(pValue2,3)	#设置p值的输出格式
#fit2=survfit(Surv(futime, fustat) ~ group2, data = rtg2)
#summary(fit1)	#查看五年生存率

#11.聚类后各组表达量差异分析#########################################################################################
nor<-as.vector(a0[1:normalNum,][,1])	#正常甲基化病人编号(只提取前四个格式便于转录组病人编号交集)
group.exp<-t(IntersectmRNAexp[gene,])	#提取该基因转录组表达矩阵

for(i in 1:3){	#第一到第三组循环
group.exp.inter<-as.vector(group.exp[intersect(rownames(group.exp),c(nor,as.vector(rt[rt$state==i,][,1]))),])	#正常对照及第i组甲基化病人的转录组表达向量
nor.inter=group.exp.inter[1:normalNum]	#提取交集正常组甲基化病人表达量向量
pat.inter=group.exp.inter[-(1:normalNum)]	#提取交集第i组甲基化病人表达量向量

factorgroup<-factor(c(rep(1,normalNum),rep(2,(length(group.exp.inter)-normalNum))))	#创建因子进行方差齐性检验
rand<-rnorm(length(group.exp.inter),sd=1e-6)	#避免数据出现组内都为0时无法进行levene.test的情况以及数据重复导致的wilcox.test"无法精确计算带连结的p值"错误,构建等宽随机数向量
group.exp.inter<-group.exp.inter+rand	#矩阵加上随机数向量
if(leveneTest(group.exp.inter~factorgroup)[[3]][1]>0.05){	#levene正态性检验,若符合正态性
fcpv<-t.test(nor.inter,pat.inter,alternative = "two.sided")[[3]]	#student's t检验
}else{	#若不符合正态性
fcpv<-wilcox.test(nor.inter,pat.inter)[[3]]	#wilcox秩和检验的p值
}
log2fc<-log(mean(pat.inter)/mean(nor.inter),2)	#计算表达量差异倍数
eval(parse(text=paste(paste("group",i,"fcpv",sep=""), '=fcpv')))	#将第一组到第三组的p值分别写入对应变量中
eval(parse(text=paste(paste("group",i,"log2fc",sep=""), '=log2fc')))	#将第一组到第三组的倍数值分别写入对应变量中
}
}
#12.数据输出#########################################################################################
Methylog2fc<-as.vector(outTab[outTab$gene==gene,4])	#查找该基因甲基化程度在正常和癌中的倍数差异
Methyfcfdr<-as.vector(outTab[outTab$gene==gene,6])	#查找该基因甲基化程度在正常和癌中倍数差异的p值

write.table(paste(gene,cor,pval,Methylog2fc,Methyfcfdr,RNAlog2fc,RNAfcpv,num1,num2,num3,pValue1,pValue2,group1log2fc,group1fcpv,group2log2fc,group2fcpv,group3log2fc,group3fcpv,sep="\t"),file="result.txt",row.names=F,col.names=F,append=TRUE,quote=F)
#library(pheatmap)
#pheatmap(a1[,-c(1,ncol(a1)-1,ncol(a1)],cluster_rows = FALSE,show_rownames=F,gaps_row=45)
}
