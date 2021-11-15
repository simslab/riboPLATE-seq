#!/usr/bin/env Rscript

suppressMessages(library("DESeq2",quietly=TRUE))
suppressMessages(library("BiocParallel",quietly=TRUE))
register(MulticoreParam(12))
args <- commandArgs(trailingOnly=TRUE)
for (i in 1:length(args))
{
        print(args[i])
}
#[1:directory 2:ID 3:outpfx 4:design 5: base1 6: base2]
inDir = args[1]
countData <- read.table(paste0(inDir,args[2],'.cts.txt'))
colData <- read.csv(paste0(inDir,args[2],'.cols.txt'),colClasses="factor")
colData$type <- relevel(colData$type,"RNA")
colData$condition <-relevel(colData$condition,args[5])
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=as.formula(args[4]))#args[4])
des = args[4]
desR = paste0(head(strsplit(des,'+',fixed=TRUE)[[1]],-1),collapse='+')
print(des)
print(desR)
if (desR == ''){
	desR <- "~1"
}
print(desR)
print(des)

for (c_base in c(args[5],args[6]))
{
	dds$condition <- relevel(dds$condition,c_base)
	for (c in c("PP_6","PP_30"))
	{
		dds_sub <- subset(dds,select=(colData(dds)$condition %in% c(c,c_base)))#==c | dds$condition==args[5],]
		dds_sub$condition <- droplevels(dds_sub$condition)
		dds_sub <- DESeq(dds_sub,parallel=TRUE,test='LRT',fitType='local', reduced=as.formula(desR))
		for (resN in resultsNames(dds_sub))
		{
			print(resN)
			res <- results(dds_sub,parallel=TRUE,name=resN)
			res <- res[order(res$padj),]
			write.table(res, paste0(inDir,args[2],'.',args[3],'.',resN,'.',c_base,'.results.txt'), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
		}
	}
}

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=as.formula(des))#~1)#args[4])
dds$condition<-relevel(dds$condition,args[5])
dds <- DESeq(dds,parallel=TRUE,fitType='local',test='LRT',reduced=as.formula(desR))#full=m1,test='Wald',fitType='local')
ncs <- counts(dds,normalized=TRUE)
write.table(assay(ncs),file=paste0(inDir,args[2],'.',args[3],'.normB.txt'),quote=FALSE,row.names=TRUE,col.names=TRUE,sep='\t')
ncs <- normTransform(dds)
write.table(assay(ncs),file=paste0(inDir,args[2],'.',args[3],'.norm.txt')
vst <- vst(dds, blind=FALSE)
write.table(assay(vst),file=paste0(inDir,args[2],'.',args[3],'.vst.txt'),quote=FALSE,row.names=TRUE,col.names=TRUE,sep='\t')