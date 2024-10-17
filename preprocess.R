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
setwd("D:/PhD/My Papers/Paper2/code/R")
#sigcpgs1 <- read.csv.ffdf(file="D:\\PhD\\methylationPaper2\\R\\Data\\GSE80417.csv",header=TRUE, VERBOSE=TRUE, first.rows=10000, next.rows=50000)
#sigcpgs1 <- read.table("D:\\PhD\\methylationPaper2\\R\\Data\\GSE144858.csv",sep=",",header=T)
sigcpgs1 <- fread("D:\\PhD\\My Papers\\Paper2\\methylationPaper2\\R\\Data\\GSE128235.txt",sep="\t",header=T)

dim(sigcpgs1)

#colnames(sigcpgs1)
mydata <- as.data.frame(sigcpgs1)
columns_to_remove <- grep("_DetectionPval", names(mydata))
mydata=mydata[,-columns_to_remove]
mydata <- mydata[,-1]
mydata2 <- mydata[,-1]
rownames(mydata2) <- mydata[,1]
#idx = which(names(mydata)=="*Detection.pval")
#Dasan Normaization
#label=as.factor(sub(".*_", "", names(mydata2)))

#names(mydata2)=as.factor(sub(".*_", "", names(mydata2)))
#levels(label) <- c("I","II")
#label=as.character(label)
#levels(names(mydata2)) <- c(as.roman(1), as.roman(2))
#mydata <- mydata[ , -574]
#mydata <- mydata[ , -691]

#Dasan normalization
#beta=dasen(as.matrix(mydata2),NULL,label)


#remove missing and duplicate values
mydata_refined<-na.omit(mydata2)
sum(duplicated(rownames(mydata_refined)))
rows=rownames(mydata_refined)

#Get phenoData
GSE111629=getGEO("GSE111629",GSEMatrix = TRUE, getGPL= FALSE)
#GSE42861 <- getGEO(filename="D:\\PhD\\methylationPaper2\\R\\Data\\GSE42861_series_matrix.txt",GSEMatrix = FALSE,getGPL = FALSE) #Retrieve matrix data and store it in R object
phenoData = as.data.frame(GSE111629[[1]]) 
#rownames(phenoData)=GSE128235[[1]]$geo_accession

#phenoData=fread("D:\\PhD\\methylationPaper2\\R\\Data\\phenoData_gse42861.csv",sep=",",header=T)
phenoData$diagnosis.ch1=as.factor(as.integer(phenoData$diagnosis.ch1 == "case"))
#phenoData$disease.state.ch1=as.factor(as.integer(phenoData$disease.state.ch1 == "Parkinson's disease (PD)"))
#phenoData[1]=as.factor(as.integer(phenoData[1] == "Smoker"))
#rownames(phenoData)=phenoData$Accession
#rownames(phenoData)=phenoData$geo_accession
#phenoData=phenoData[match(colnames(mydata_refined),paste0("X",phenoData$description)),]
#phenoData=phenoData[match(colnames(mydata_refined),rownames(phenoData)),,drop=FALSE]
#phenoData=phenoData[match(colnames(mydata_refined),phenoData$source_name_ch1),]
phenoData$title=paste0("sample",sub(".* ", "", phenoData$title))
phenoData=phenoData[match(colnames(mydata_refined),phenoData$title),,drop=FALSE]
phenoData=phenoData[rowSums(is.na(phenoData)) != ncol(phenoData), ]
colnames(mydata_refined)=rownames(phenoData)


#get probGene mapping
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_data<-ann[c("Name","UCSC_RefGene_Name")]
ann_DF <-as.data.frame(ann_data)
ann_DF<-ann_DF[!(is.na(ann_DF$UCSC_RefGene_Name) | ann_DF$UCSC_RefGene_Name==""), ]
rownames(ann_DF)<-ann_DF[,1]
ann_DF_sub<-ann_DF[match(rownames(mydata_refined),rownames(ann_DF)),]
ann_DF_sub=ann_DF_sub[rowSums(is.na(ann_DF_sub)) != ncol(ann_DF_sub), ]
ann_DF_sub<-ann_DF_sub[!(is.na(ann_DF_sub$UCSC_RefGene_Name) | ann_DF_sub$UCSC_RefGene_Name==""), ]
ann_DT<-setDT(ann_DF_sub)
df1<-cSplit(ann_DT, c("UCSC_RefGene_Name"), sep = ";", direction = "long")
dim(df1)
df2<-na.omit(df1, cols=c("UCSC_RefGene_Name"))
dim(df2)
df3<-unique(df2, by = c('Name', 'UCSC_RefGene_Name'))
dim(df3)
mydata_refined_sub<-mydata_refined[match(unique(df3$Name),rownames(mydata_refined)),]
mydata_refined_sub<-na.omit(mydata_refined_sub)
dim(mydata_refined_sub)


############################################
fwrite(phenoData,"phenoData_GSE128235.csv",row.names=TRUE)
fwrite(mydata_refined_sub,"processedData_GSE128235.csv",row.names=TRUE)
fwrite(df3, "probeGene_GSE128235.csv",row.names=TRUE)


