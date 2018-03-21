#Rscript --vanilla $Heatmap $folder $seg $min $max $sample_annotation
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(data.table)
library(readr)
library(reshape2)
library(grid)

args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("folder","seg","min","max","sample_annotation")
print(args)

segment=args["seg"]
folder=args["folder"]
annotation <- read.table(args["sample_annotation"], header=TRUE)
st=seq(as.numeric(args["min"]),as.numeric(args["max"]), as.integer( ( as.numeric(args["max"]) - as.numeric(args["min"]) )/( as.numeric(args["min"]) -1 ) ) )


#read MethylSeekR output and prepare the data to plot average methylation
temp <- list.files(paste(folder,segment,"s/bed/",sep=""), full.names = TRUE)
files = lapply(temp,function(x){read.table(gzfile(x),skip=1)})
d1 <- as.data.frame(rbindlist(files))
colnames(d1) <- c("chr","start","end","segments","meth","direction","start","end","color")
names <- lapply(temp, function(x) {str <- basename(x); spl <- strsplit(str,"\\."); paste(spl[[1]][1:1],collapse="_")} )
names.vector <- c()
for(i in 1:length(names))
{
names.vector <- c(names.vector,rep(names[[i]][1],length(files[[i]]$V1)))
}

d1$sample <- names.vector
d1$length <- d1$end-d1$start

#prepare the data for stacked bar plots
tmp <- data.frame(d1[,c(4,10,11)])
m <- melt(tmp,id.vars=c("segments","sample"))
d <- aggregate(formula=value ~ sample + segments , data=m, FUN= sum)
d2 <- dcast(d,sample ~ segments)
d2[,2:5] <-d2[,2:5]/rowSums(d2[2:5])
d3 <-melt(d2)
colnames(d3) <- c("sample","segments","length")
d3[,"segments"] <- factor(d3[,"segments"], levels=c("FMR","PMD","LMR","UMR"))

#plot the heatmap and reordered box- and stackedBar- plots
for(i in st){
        #Heatmap
        cat("plotting Heatmap for",i,"states of",segment,"\n")
        file=paste(folder,segment,"s/",segment,"s_chromhmm/",segment,"s_mix_",i,"St/emissions_",i,".txt",sep="")
        print(file)
        x <- read.delim(file)
        name=paste(folder,segment,"s/",segment,"s_chromhmm/",segment,"s_mix_",i,"St/emission_",i,"_",segment,"_","HC.svg",sep="")
        print(name)
        df <- x[,-1]
        rownames(df) <- x[,1]
        names(df) <- gsub("^X","",names(df))
        colnames(df)[colnames(df)=="Plas_BM_bmPCs.V156"] <- "Plas_BM_bmPCs_V156"
        annotationOrd <- merge(data.frame(sample=colnames(df)), annotation, sort=FALSE)
        ann <- HeatmapAnnotation(annotationOrd[,-c(1)])

        hm <- Heatmap(df, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), top_annotation=ann,clustering_method_columns="ward.D2",column_names_gp=gpar(cex=0.5), name= segment, column_title = segment, column_dend_height = unit(4, "cm"))
        svg(name, width=17, height=9)
        print(hm)
        dev.off()
        
        #boxplots
        cat("plotting boxplots for",i,"states of",segment,"\n")
        d1.levels <- names(df)[column_order(hm)]
        color_order <- ann@anno_list[[colnames(annotation)[3]]]@color_mapping@colors
        d1_new <- merge(d1[,c("segments","meth","sample","length")], annotationOrd, sort=FALSE)
        d1_new[,'sample'] <- factor(d1_new[,'sample'], levels=d1.levels)
        name_boxplot = paste(folder,segment,"s/",segment,"s_chromhmm/",segment,"s_mix_",i,"St/avg_meth_",i,"_",segment,".png",sep="")
        print(name_boxplot)

        gg1 <- ggplot(subset(d1_new,segments==segment), aes_string(x="sample", y="meth", fill=colnames(annotation)[3])) + geom_boxplot(outlier.size=0.3)  + theme(legend.position = "none", strip.text.x = element_text(colour = "white"),strip.background = element_rect(colour = "white", fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Arial Black"), axis.text.x=element_text(family="Arial Black",colour = "black",  size=6, angle=90), axis.ticks.x=element_blank()) + xlab("") + ylab("weighted average methylation") + scale_fill_manual(values=color_order)

        ggsave(filename=name_boxplot, plot=gg1, device="png", width=17, height=9)
        
        #stacked barplots
        cat("plotting stacked barplots for",i,"states of",segment,"\n")
        d3_new <- merge(d3, annotationOrd, sort=FALSE)
        d3_new[,'sample'] <- factor(d3_new[,'sample'], levels=d1.levels)
        name_bar = paste(folder,segment,"s/",segment,"s_chromhmm/",segment,"s_mix_",i,"St/segments_perc_",i,"_",segment,".png",sep="")
        print(name_bar)

        gg2 <- ggplot(d3_new, aes(x=sample,y=length, fill=segments)) + geom_bar(stat="identity") + theme(legend.position = "top", strip.text.x = element_text(colour = "white"),strip.background = element_rect(colour = "white", fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Arial Black"), axis.text.x=element_text(family="Arial Black",colour = "black",  size=6, angle=90)) + ylab("genome percentage %") + scale_fill_manual(values=c("#CA0020","#F4A582","#92C5DE","#0571B0"))

        ggsave(filename=name_bar, plot=gg2, device="png", width=17, height=9)
}

