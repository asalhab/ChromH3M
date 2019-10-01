
set.seed(123)
library(MethylSeekR)
library(ggplot2)

args <- commandArgs(trailingOnly= TRUE)
names(args) <- c("genome","WD_output","file_Name","snp_file","Sample_Name","prefix","meth.level","FDR","minCover","nrCores")
print(args)

setwd(args["WD_output"])

if (args["genome"]=="hg19") {
cat("extracting CpG islands from UCSC ......\n")
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg19")
sLengths=seqlengths(Hsapiens)
myGenomeSeq <- Hsapiens
session <- browserSession(url = 'http://genome-euro.ucsc.edu/cgi-bin/')
cat("loading session ......\n")
genome(session) <- "hg19"
cat("make a query ......\n")
query <- ucscTableQuery(session, "cpgIslandExt")
cat("make the track .......\n")
#CpGislands.gr <- track(query, asRangedData = FALSE)
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
cat("supressWarnings ......\n")
CpGislands.gr <-suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))
cat("CpG islands have been extracted ......\n")
}

if(args["genome"]=="hg38"){
cat("extracting CpG islands from UCSC ......\n")
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
sLengths=seqlengths(Hsapiens)
myGenomeSeq <- Hsapiens
session <- browserSession(url = 'http://genome-euro.ucsc.edu/cgi-bin/')
cat("loading session ......\n")
genome(session) <- "hg38"
cat("make a query ......\n")
query <- ucscTableQuery(session, "cpgIslandExt")
cat("make the track .......\n")
#CpGislands.gr <- track(query, asRangedData = FALSE)
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
cat("supressWarnings ......\n")
CpGislands.gr <-suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))
cat("CpG islands have been extracted ......\n")

}

if(args["genome"]=="mm10") {
library(rtracklayer)
library("BSgenome.Mmusculus.UCSC.mm10")
sLengths=seqlengths(Mmusculus)
myGenomeSeq <- Mmusculus
session <- browserSession(url = 'http://genome-euro.ucsc.edu/cgi-bin/') ## a url was added here becasue browserSession() didn't work for mm10
#session <- browserSession()
cat("loading session ......\n")
genome(session) <- "mm10"
cat("make a query ......\n")
query <- ucscTableQuery(session, "cpgIslandExt")
cat("make the track .......\n")
#CpGislands.gr <- track(query, asRangedData = FALSE)
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
cat("supressWarnings ......\n")
CpGislands.gr <-suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))
cat("CpG islands have been extracted ......\n")
}

### read input file 
cat("Reading the input file ......\n")
data.gr <- readMethylome(args["file_Name"], seqLengths=sLengths)

### read SNPs
cat("Reading the snp file ......\n")
snps.gr <- readSNPTable(args["snp_file"], seqLengths=sLengths)

### removing SNP
cat("Filtering the SNPs ......\n")
data.gr <- removeSNPs(data.gr, snps.gr)

cat("The mean coverage of this sample after removing SNPs is:", mean(data.gr$T), "\n")
cat("The median coverage of this sample after removing SNPs is:", median(data.gr$T), "\n")
### plot alpha distribution
cat("Plotting alpha distribution ......\n")
plotAlphaDistributionOneChr(m=data.gr, chr.sel="chr2",pdfFilename=paste(args["Sample_Name"],args["prefix"],"alpha.distribution.pdf",sep="."),num.cores=as.numeric(args["nrCores"]))

### find PMD and plot
cat("Finding PMDs ......\n")
PMDsegments.gr <- segmentPMDs(m=data.gr,pdfFilename=paste(args["Sample_Name"],args["prefix"],"alpha.model.fit.pdf",sep="."), chr.sel="chr2",seqLengths=sLengths, num.cores=as.numeric(args["nrCores"]))
plotAlphaDistributionOneChr(m=subsetByOverlaps(data.gr,PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="chr2",pdfFilename=paste(args["Sample_Name"],args["prefix"],"alpha.distribution.without.PMDs.pdf",sep="."),num.cores=as.numeric(args["nrCores"]))
savePMDSegments(PMDs=PMDsegments.gr,GRangesFilename=paste(args["Sample_Name"],args["prefix"],"PMDs.gr.rds",sep="."), TableFilename=paste(args["Sample_Name"],args["prefix"],"PMDs.tab",sep="."))

### FDR calculation
cat("FDR calculation ......\n")
stats <- calculateFDRs(m=data.gr, CGIs=CpGislands.gr,PMDs=PMDsegments.gr, nCpG.cutoffs =seq(1, 17, by=3),pdfFilename=paste(args["Sample_Name"],args["prefix"],"FDR.pdf",sep="."),num.cores=as.numeric(args["nrCores"]))
FDR.cutoff <- as.numeric(args["FDR"])
m.sel <- as.numeric(args["meth.level"])
n.sel=as.integer(names(stats$FDRs[as.character(m.sel), ][stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
cat("minimum number of CpGs in LMRs:",n.sel,"CpGs\n")
saveRDS(stats, file=paste(args["Sample_Name"],args["prefix"],"stats.rds",sep=".") )

### find UMR and LMR
cat("Finding UMRs and LMRs ......\n")
UMRLMRsegments.gr <- segmentUMRsLMRs(m=data.gr, meth.cutoff=m.sel,nCpG.cutoff=n.sel, PMDs=PMDsegments.gr,pdfFilename=paste(args["Sample_Name"],args["prefix"],"UMR.LMR.scatter.plot.pdf",sep="."),num.cores=as.numeric(args["nrCores"]), myGenomeSeq=myGenomeSeq,seqLengths=sLengths, minCover = as.numeric(args["minCover"]))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=paste(args["Sample_Name"],args["prefix"],"UMRsLMRs.gr.rds",sep="."),TableFilename=paste(args["Sample_Name"],args["prefix"],"UMRsLMRs.tab",sep="."))

### plot the final segmentation
cat("Plotting the final segmentation ......\n")
plotFinalSegmentation(m=data.gr, segs=UMRLMRsegments.gr,PMDs=PMDsegments.gr,numRegions = 4,pdfFilename=paste(args["Sample_Name"],args["prefix"],"final.segmentation.example.regions.pdf",sep="."),meth.cutoff=m.sel)

### plot the methylation and coverage distibutions
df <- as.data.frame(data.gr)
df$meth <- df$M/df$T
pdf(paste(args["Sample_Name"],args["prefix"],"meth.cov.pdf",sep="."))
ggplot(df[df$T>=5,], aes(x=df[df$T>=5,"meth"])) + geom_density(colour="dodgerblue1",size=1) +  ylab("density") + xlab("beta value") + ggtitle("Methylation level density")
ggplot(df, aes(x=df$T)) + geom_histogram(binwidth=1,alpha=.5,position="identity",colour = "dodgerblue1", fill = "dodgerblue1") + geom_vline(xintercept=mean(df$T),colour="black", linetype = "longdash")  + ylab("Frequency") + xlab("read coverage per CpG") + geom_text(aes(x2,y2,label = texthere),data.frame(x2=mean(df$T), y2=max(table(df$T)), texthere=round(mean(df$T),2)))
dev.off()
