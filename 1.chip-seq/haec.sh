######################## Variate ################################
dir="/home/lqd/GENCODE/"
######################## Tools ##################################
fastqc="/root/QualityControl/FastQC/fastqc"
sratools="/root/sratoolkit.2.7.0-ubuntu64/bin/fastq-dump"
filter="/root/QualityControl/fastx_toolkit-0.0.14/src/fastq_quality_filter/fastq_quality_filter -q 20 -p 70"
bbduk="/home/lqd/Software/bbmap/bbduk.sh"
repair="/home/lqd/Software/bbmap/repair.sh"
ref="/home/lqd/Software/bbmap/resources/adapters.fa"
bowtie="/root/Mapping/bowtie2-2.2.9/bowtie2 --met-file ./record.txt -p 15 -x /home/lqd/Database/hg38/hg38"
macsparameters="-g hs -p 1e-5 -B -S"
ToBigWig="/home/lqd/Software/ucsctoolkits/bedGraphToBigWig"
chromsizes="/home/lqd/Software/ucsctoolkits/hg38.chrom.sizes"
######################## Tools ##################################
for cell in HAEC
  do
    tag=""$cell"_"
	inputdir="$dir$cell/"
	#mkdir "$inputdir"fastqc
    for treat in H3K27ac
      do
	filename="$tag$treat"
	cd $inputdir
	#$fastqc $filename*.fastq -t 14 -o ./fastqc
        #$filter -i $filename*.fastq -o $filename.clean1.fastq
        #$fastqc $filename.clean1.fastq -t 14 -o ./fastqc
      done
  done
  for cell in HAEC
  do
    tag=""$cell"_"
	inputdir="$dir$cell/"
    for treat in H3K27ac
      do
	#filename="$tag$treat"
	#cd $inputdir
        #$bowtie $filename.clean1.fastq -S $filename.sam	
	#grep -v "XS:i:" $filename.sam | samtools view -bS - > $filename.unique.bam
      done
  done
for cell in HAEC
  do
    tag=""$cell"_"
    for treat in H3K27ac
      do
	inputdir="$dir$cell/"
        filename="$tag$treat"
	cd $inputdir
	#macs14 -t $filename.unique.bam $macsparameters -n $treat
	cd "$inputdir$treat"_MACS_bedGraph/treat/
	gunzip *.gz
	sum=$(samtools view "$inputdir$filename".unique.bam | awk 'BEGIN{sum=0}{if($2==16 || $2==0){sum +=1}};END{print sum}')
	cat *treat_afterfiting*.bdg | awk '{print $1"\t"$2"\t"$3"\t"$4*1000000/"'$sum'"}' > "$cell"_$treat.bdg
	$ToBigWig "$cell"_$treat.bdg $chromsizes "$cell"_$treat.bw
	rm "$cell"_$treat.bdg
      done
  done