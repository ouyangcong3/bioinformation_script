main="/home/lqd/GENCODE/HCT116/ColonGeneExpression/TCGA/DNAmethylation/" #该文件夹下需要.json文件,survival_time.pl,genename.txt,normalizeExp.txt,TSS2kbProbespos.txt
dir="/home/lqd/GENCODE/HCT116/ColonGeneExpression/TCGA/DNAmethylation/SourceFile/" #该文件下存放下载好的TCGA甲基化数据

#1.Move the '.txt' file from each folder
cd $dir
for file in ` ls $1 `
  do
    cd $file
    mv *.txt $dir
    cd $dir
  done

#2.Extract the rank of patients id from the '.josn' file  
cd $main
cat *.json | awk '{if($1 ~/^"file_name/){print $0}}' > filename.txt
cat filename.txt | awk -F '.' '{print $6}' > TCGAid.txt

#3.Sort the rank by patten 'Tumor-Control'
cat TCGAid.txt | awk -F '-' '{if($4 ~/^11/){print $0}}' > Control.txt
cat TCGAid.txt | awk -F '-' '{if($4 ~/^11/){}else{print $0}}' > Tumor.txt
cat Control.txt Tumor.txt  > TCGAfilename.txt
cat TCGAfilename.txt | awk -F '-' '{if($4 ~/^01/){print $0"\t""T"}}{if($4 ~/^11/){print $0"\t""C"}}' > group.txt

#4.Combind the Matrix
file=$(cat TCGAfilename.txt | head -1) #提取列表中第一个病人的病人编号
cat $dir*$file*.txt | tail -n +2 | awk '{print $3"\t"$4"\t"$5"\t"$1}' > 450Kprobes.bed #制作450K甲基化微阵列的探针集,提取第3-5列和第1列制作为bed格式文件
cat $dir*$file*.txt | tail -n +2 | awk '{print $6}' > probesann.txt #提取甲基化注释信息
cat $dir*$file*.txt | awk '{print $1"\t"$2}' > temporary.txt #读取上一步所得病人编号的文件,提取前两列(探针号和beta值)
sed '1c id\t'$file'' temporary.txt > Initiation.txt #将第一行修改为"id\t病人编号"格式并输出到Initiation.txt文件
for file in ` cat TCGAfilename.txt | tail -n +2 ` #循环提取列表中其它病人的编号
  do
    cat $dir*$file*.txt | awk '{print $2}' > temporary.txt #读取该编号的文件,提取第二列(beta值)
    sed '1c '$file'' temporary.txt > $file #将第一行改为病人编号并输出到病人编号命名的变量
  done
paste Initiation.txt ` tail -n +2 TCGAfilename.txt ` > MethylationMatrix.txt #按TCGAfilename.txt中的顺序合并表格并输出
rm TCGA-* #删除临时文件
rm temporary.txt Initiation.txt Control.txt Tumor.txt TCGAid.txt filename.txt #删除临时文件

##5.Check the Matrix
#head -1 MethylationMatrix.txt | awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' > 1.txt
#head -2 MethylationMatrix.txt | awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' > 2.txt

##6.Make TSS2kbProbes file (First download TSS2kb.hg38.bed from UCSC)
#bedtools intersect -a 450Kprobes.bed -b RefSeqTSS2kb.hg38.bed > TSS2kbProbespos.txt #将450K甲基化微阵列的探针集中位于TSS2kb范围内的探针提取出来
#cat TSS2kbProbespos.txt | awk '{print $4}' > TSS2kbProbes.txt #提取表中的探针号

#7.Prepare survival time
perl survival_time.pl *.json

#8.Transform
Rscript Methy-Transcriptome-Analysis.r
