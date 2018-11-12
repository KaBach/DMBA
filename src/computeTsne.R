# Estimate size factors using scran
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(ggplot2) 
source("../../../LinTracing/repo/src/functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
	      
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

# Gene and cell filtering
m <- m[,pD$PassAll]
fD$keep <- rowMeans(m)>0.01 
pD <- pD[pD$PassAll,]
m <- m[fD$keep,]
fD <- fD[fD$keep,]
rownames(m) <- as.vector(rownames(m))

clusters <- quickCluster(m,method="igraph")
minSize <- min(table(clusters))
param <- MulticoreParam(workers=4)
pD$sf <- computeSumFactors(m, sizes=seq(20,min(100,minSize),5),clusters=clusters,positive=TRUE, BPPARAM=param)

plot(log10(colSums(m))~log10(pD$sf),main="Library Size versus Size Factors")
#
out <- pD[,c("barcode","sf")]
write.table(file="../data/Robjects/SizeFactors.csv", out, row.names=FALSE, sep=",")

m <- log2(t(t(m)/pD$sf)+1)

# Highly variable genes
hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
hvg <- hvg[rownames(hvg) %in% fD$id[fD$KeepForHvg],]
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

# Compute tSNE 
fPCA <- m[hvg,]
fPCA <- sqrt((1-cor(as.matrix(fPCA)))/2)
set.seed(300)
tsn <- Rtsne(fPCA,perplexity=30,check_duplicates=FALSE,is_distance=TRUE)
pD$tSNE1 <- tsn$Y[,1]
pD$tSNE2 <- tsn$Y[,2]
write.csv(file="../data/Robjects/tSNE.csv",pD[,c("barcode","tSNE1","tSNE2")])
