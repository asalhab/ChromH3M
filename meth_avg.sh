#!/bin/bash
#$ -S /bin/bash
#$ -N meth_avg
#$ -cwd
#$ -j y
#$ -V
#$ -l h_vmem=40G

Pipeline_name="multi_Seg_summary"

printHelp() {
   echo -e "$(tput bold)Description:$(tput sgr0)"
   echo -e ""
   echo -e "$(tput bold)Usage:$(tput sgr0)"
   echo -e " bash $Pipeline_name  $(tput bold)$(tput setaf 1)-i files_dir    -s segment     -o output_dir     -n sample name     -g genome$(tput sgr0)"
   echo -e ""
   echo -e " $(tput bold)Mandatory:$(tput sgr0)"
   echo -e "  -i directory contains gzipped segments files"
   echo -e "  -s segment (PMD, FMR ...)"
   echo -e "  -o output folder"
   echo -e "  -n Sample name"
   echo -e "  -g genome length (Shortcuts: hg19, mm10, mm9)"
}

while getopts ":hi:s:o:n:g:" opt
do
   case $opt in
      h) printHelp; exit 0 ;;
      i) files_dir="$OPTARG" ;;
      s) segment="$OPTARG" ;;
      o) out_dir="$OPTARG" ;;
      n) name="$OPTARG" ;;
      g) genome="$OPTARG" ;;
      *) printHelp; exit 1 ;;
   esac
done

R="directory to R"

if [[ -z "$files_dir" ]] || [[ -z "$out_dir" ]] || [[ -z "$name" ]] || [[ -z "$genome" ]]
   then
   echo -e "\n$(tput bold)$(tput setaf 1)ERROR ($Pipeline_name): Must set all mandatory options$(tput sgr0)\n"
   printHelp
   exit 1
fi


case "$genome" in
    hg19) genome="Path to the human genome length file (chr, length)"
          assembly="hg19"
          gaps="Path to the bed file containing the gaps you want to filter (gaps from UCSC)"
          ;;

    mm10) genome="Path to the mouse genome length file (chr, length"
          assembly="mm10"
          gaps="Path to the bed file containing the gaps you want to filter (gaps from UCSC"
         ;;

    mm9) genome="Path to the mouse genome length file (chr, length)"
         assembly="mm9"
         gaps="Path to the bed file containing the gaps you want to filter (gaps from UCSC)"
         ;;
esac

case "$segment" in
    PMD) filter=10000
         window=1000
         step=1000
         ;;
    FMR) filter=2000
         window=1000
         step=1000
         ;;
    LMR) filter=200
         window=200
         step=200
         ;;
    UMR) filter=200
         window=200
         step=200
         ;;
esac


###### generating files contating the overlap of each sample's PMDs(>10kb) with 1kb binned genome and then merging them to be analysed by clustering, showing the similarity between PMDs across samples
echo -e "generating the ${segment} files ......."

for i in $files_dir/*gz
do
        zcat $i|awk 'NR>1'|sed 's/chr//'|grep $segment > $out_dir/$(basename ${i%bed.gz}${segment}.bed)
done



echo -e "generating the binned genome file ......"
bedtools makewindows -g <(egrep -v "GL|JH|phi|L|M|hs|NC|Y" $genome) -w $window -s $step > $out_dir/${assembly}.windows$((window/1000))kb.bed

echo -e "generating the counts files ......"
for i in $out_dir/*${segment}.bed
do
        bedtools intersect -a $out_dir/${assembly}.windows$((window/1000))kb.bed -b <(awk -v f=$filter '($3-$2)>f' $i|grep -v "NA") -c \
        | awk -vOFS='\t' '{print $1,$2,$3,$4}' \
        | bedtools merge -i - -c 4 -o distinct -d -1 \
        > $out_dir/$(basename ${i%bed}counts.bed); done


cd $out_dir
names=(*${segment}.counts.bed)


paste $(for i in `seq 0 $(( ${#names[@]} -1 ))`; do echo ${names[$i]};done|tr '\n' ' ') |cut -f 1-3,$(seq  4 4 $(( ${#names[@]} * 4 )) |tr '\n' ',' |sed s/,$//)  - |cat <(echo "chr start end ${names[@]%%.*}"| tr ' ' '\t') -  > merged.sampels.${segment}.bed


echo -e "starting R ......"

### R commands
$R  --vanilla <<RSCRIPT
library(ggplot2)
library(data.table)
library(readr)
library(reshape2)

setwd("$out_dir")
counts <- read.table("$out_dir/merged.sampels.${segment}.bed",header=TRUE)
df <- as.data.frame(counts[counts\$chr!="X",4:length(counts)])
hc <- hclust(dist(t(df)),"complete")
png(paste("$name",".${segment}.similarity.complete.png",sep=""),width=800,height=600)
plot(hc,hang=-1,main=paste("$name","${segment} similarity",sep=" "))
dev.off()

hc <- hclust(dist(t(df)),"ward.D2")
png(paste("$name",".${segment}.similarity.wardD2.png",sep=""),width=800,height=600)
plot(hc,hang=-1,main=paste("$name","${segment} similarity",sep=" "))
dev.off()

png(paste("$name",".${segment}.similarity.histogram.png",sep=""))
hist(rowSums(counts[,-c(1,2,3)]),breaks=50)
dev.off()


temp <- list.files("$files_dir",pattern="*segments.bed.gz$",recursive=TRUE,full.names=TRUE)
files = lapply(temp,function(x){read.table(gzfile(x),skip=1)})
#df <- c()
#for(i in 1:length(files))
#{
#df <- data.frame(rbind(df,files[[i]]))
#}

df <- as.data.frame(rbindlist(files))
colnames(df) <- c("chr","start","end","segment","meth","direction","start","end","color")
names <- lapply(temp, function(x) {str <- basename(x); spl <- strsplit(str,"\\\."); paste(spl[[1]][1:1],collapse="_")} )
signals <- lapply(temp, function(x) {str <- basename(x); spl <- strsplit(str,"_"); paste(spl[[1]][5],collapse="_")} )
replicates <- lapply(temp, function(x) {str <- basename(x); spl <- strsplit(str,"_"); paste(spl[[1]][2],collapse="_")} )

names.vector <- c()
for(i in 1:length(names))
{
names.vector <- c(names.vector,rep(names[[i]][1],length(files[[i]]\$V1)))
}

signals.vector <- c()
for(i in 1:length(names))
{
signals.vector <- c(signals.vector,rep(signals[[i]][1],length(files[[i]]\$V1)))
}

replicates.vector <- c()
for(i in 1:length(names))
{
replicates.vector <- c(replicates.vector,rep(replicates[[i]][1],length(files[[i]]\$V1)))
}

df\$sample <- names.vector
df\$signal <- signals.vector
df\$replicate <- replicates.vector

df[,"segment"] <- factor(df[,"segment"], levels=c("FMR","PMD","LMR","UMR"))
df <- df[!grepl("X|Y",df\$chr),]
df <- df[complete.cases(df),]
df\$length <- df\$end-df\$start

png(paste("$name",".avg.meth.png",sep=""),width=800,height=600)
ggplot(df, aes(x=sample, y=meth, fill=segment)) + geom_boxplot(outlier.size=1)  + theme(legend.position = "top", strip.text.x = element_text(colour = "white"),strip.background = element_rect(colour = "white", fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),text=element_text(family="Arial Black"), axis.text.x=element_text(family="Arial Black",colour = "black",  size=10, angle=90)) + scale_fill_manual(values=c("#CA0020","#F4A582","#92C5DE","#0571B0")) +xlab("") +ylab("weighted average methylation") 
dev.off()

png(paste("$name",".","$segment",".avg.meth.png",sep=""),width=800,height=600)
ggplot(df[grepl("$segment",df\$segment),], aes(x=sample, y=meth, fill=segment)) + geom_boxplot(outlier.size=1)  + theme(legend.position = "top", strip.text.x = element_text(colour = "white"),strip.background = element_rect(colour = "white", fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), text=element_text(family="Arial Black"), axis.text.x=element_text(family="Arial Black",colour = "black",  size=10, angle=90)) + scale_fill_manual(values=c("#F4A582")) +xlab("") +ylab("weighted average methylation") 
dev.off()

#if(length(names)<50)
#{
#png(paste("$name",".density.meth.png",sep=""),width=800,height=600)
#ggplot(df, aes(x=length, y=meth)) + geom_point(size=0.3) + geom_density2d(aes(colour=segment)) + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_colour_manual(values=c("#CA0020","#F4A582","#92C5DE","#0571B0")) + facet_wrap(~sample) + theme(legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), text=element_text(family="Arial Black"))+ annotation_logticks(sides = "b") + xlab("length(bp)")  +ylab("average methylation")
#dev.off()
#}

tmp <- data.frame(df[,c("segment","sample","length")])
m <- melt(tmp,id.vars=c("segment","sample"))
d <- aggregate(formula=value ~ sample + segment , data=m, FUN= sum)
d2 <- dcast(d,sample ~ segment)
d2[,2:length(d2)]<-d2[,2:length(d2)]/rowSums(d2[2:length(d2)])
d3<-melt(d2)
colnames(d3) <- c("sample","segment","length")

png(paste("$name",".segment.percentage.png",sep=""),width=800,height=600)
ggplot(d3) + aes(x=sample,y=length, fill=segment) + geom_bar(stat="identity") + theme(legend.position = "top", strip.text.x = element_text(colour = "white"),strip.background = element_rect(colour = "white", fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),text=element_text(family="Arial Black"), axis.text.x=element_text(family="Arial Black",colour = "black",  size=10, angle=90)) + scale_fill_manual(values=c("#CA0020","#F4A582","#92C5DE","#0571B0")) + ylab("genome percentage %")
dev.off()

RSCRIPT

rm $out_dir/*{segments.${segment}.bed,counts.bed}
#rm $out_dir/windows$((window/1000))kb.bed

