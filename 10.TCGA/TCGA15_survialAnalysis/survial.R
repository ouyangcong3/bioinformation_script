######Video source: https://shop119322454.taobao.com
#install.packages("survival")

setwd("C:\\Users\\YunGu\\Desktop\\TCGA15_survialAnalysis")   #工作目录（需修改）
library(survival)
rt=read.table("survivalInput.txt",header=T,sep="\t")
rt$futime=rt$futime/365       #如果以月为单位，除以30；以年为单位，除以365
a=rt[,"expression"]<median(rt[,"expression"])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survival.pdf")
plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
     main=paste("surival curve (p=", pValue ,")",sep=""))
legend("topright", c("BRCA1 high expression", "BRCA1 low expression"), lty = 2:3, col=c("red","blue"))
dev.off()
