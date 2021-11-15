#!/usr/bin/env Rscript
suppressMessages(library("DESeq2",quietly=TRUE))
suppressMessages(library("BiocParallel",quietly=TRUE))
register(MulticoreParam(12))
args <- commandArgs(trailingOnly=TRUE)
for (i in 1:length(args))
{
        print(args[i])
}
#[1:directory 2:ID 3:outpfx 4:design]
inDir = args[1]
countData <- read.table(paste0(inDir,args[2],'.cts.txt'))
colData <- read.csv(paste0(inDir,args[2],'.cols.txt'),colClasses="factor")
colData$type <- relevel(colData$type,"RNA")
des <- args[4]
desR <- paste0(head(strsplit(args[4],'+',fixed=TRUE)[[1]],-1),collapse='+')
if(desR == ''){
	desR <- '~1'
}
print(des)
print(desR)
if ("condition" %in% colnames(colData))
{
	colData$condition <-relevel(colData$condition,"CTRL")
}
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=as.formula(des))#args[4])
for (c in levels(dds$condition))
{
	if (c != 'CTRL'){
		dds_sub <- subset(dds,select=((colData(dds)$condition==c) | (colData(dds)$condition=='CTRL')))
		dds_sub$condition <- droplevels(dds_sub$condition)
		print(c)
		dds_sub <- DESeq(dds_sub,parallel=TRUE,fitType='local',test='LRT',reduced=as.formula(desR))
		for (resN in resultsNames(dds_sub))
		{
			res <- results(dds_sub,parallel=TRUE,name=resN)
			res <- res[order(res$padj),]
			write.table(res, paste0(inDir, args[2], '.', args[3], '.', resN, '.', c, '.results.txt'), sep='\t', row.names=TRUE, col.names=TRUE,quote=FALSE)
		}
	}
}

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=as.formula(des))#args[4])
dds <- DESeq(dds,parallel=TRUE,fitType='local',test='LRT',reduced=as.formula(desR))
ncs <- normTransform(dds)
write.table(assay(ncs),file=paste0(inDir,args[2],'.',args[3],'.norm.txt'),quote=FALSE,row.names=TRUE,col.names=TRUE,sep='\t')
vst <- vst(dds, blind=FALSE)
write.table(assay(vst),file=paste0(inDir,args[2],'.',args[3],'.vst.txt'),quote=FALSE,row.names=TRUE,col.names=TRUE,sep='\t')