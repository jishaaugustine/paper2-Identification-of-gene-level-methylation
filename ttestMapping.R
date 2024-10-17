## Not run:  # to avoid timeout on Bioconductor build
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(org.Hs.eg.db)
library(limma)
#library(missMethyl)
library(dplyr)
library(ff)
library(ffbase)
library(splitstackshape)
library(GEOquery)
library(data.table)
library(wateRmelon)
setwd("D:/PhD/methylationPaper2/R")
sUMAP <- fread("D:\\PhD\\methylationPaper2\\R\\GSE128235sUMAP.txt",sep="\t",header=T)
site <- fread("D:\\PhD\\methylationPaper2\\R\\GSE128235site.txt",sep="\t",header=T)

dim(sUMAP)
dim(site)

#get probGene mapping
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_data<-ann[c("Name","UCSC_RefGene_Name")]
ann_DF <-as.data.frame(ann_data)
ann_DF<-ann_DF[!(is.na(ann_DF$UCSC_RefGene_Name) | ann_DF$UCSC_RefGene_Name==""), ]
rownames(ann_DF)<-ann_DF[,1]
ann_DF_sub<-ann_DF[match(site$residue,rownames(ann_DF)),]
ann_DF_sub=ann_DF_sub[rowSums(is.na(ann_DF_sub)) != ncol(ann_DF_sub), ]
ann_DF_sub<-ann_DF_sub[!(is.na(ann_DF_sub$UCSC_RefGene_Name) | ann_DF_sub$UCSC_RefGene_Name==""), ]
ann_DT<-setDT(ann_DF_sub)
df1<-cSplit(ann_DT, c("UCSC_RefGene_Name"), sep = ";", direction = "long")
dim(df1)
df2<-na.omit(df1, cols=c("UCSC_RefGene_Name"))
dim(df2)
df3<-unique(df2, by = c('Name', 'UCSC_RefGene_Name'))
dim(df3)
site_sub<-merge(x = site, y = df3, by.x = "residue",by.y="Name", all.x = TRUE)
site_sub$pvalue=-log10(site_sub$pvalue)
sUMAP_sub=sUMAP[1:100,]
sUMAP_sub$pvalue=-log10(sUMAP_sub$pvalue)
d=data.frame()
for (ind in seq(0:(length(sUMAP_sub$gene)-1)))
{
  gene=sUMAP_sub[ind]$gene
  p1=sUMAP_sub[ind]$pvalue
  site_sub1=site_sub[site_sub$UCSC_RefGene_Name== gene,]
  site_sub2=site_sub1[which.max(site_sub1$pvalue),]
  d1=c(gene,p1,site_sub2$pvalue)
  d<-rbind(d,d1)
}
colnames(d)=c("gene","p1","p2")


fwrite(d,"ttest_GSE128235.csv",row.names=TRUE)


################################################################
library(ggplot2)
library(reshape)
Molten <- melt(d, id.vars = "gene")
Molten$value=as.numeric(Molten$value)
p<-ggplot(Molten, aes(x = gene, y = value, colour = variable))+ geom_point()
p + scale_y_continuous(limits=c(2,8), breaks=c(2, 4, 6, 8))


library(ggplot2)
library(reshape)
data <- data.frame(time = seq(0, 23), noob = rnorm(24), plus = runif(24), extra = rpois(24, lambda = 1))
Molten <- melt(data, id.vars = "time")
ggplot(Molten, aes(x = time, y = value, colour = variable)) + geom_line()
