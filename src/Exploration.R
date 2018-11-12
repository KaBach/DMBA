# Estimate size factors using scran
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(ggplot2) 
library(viridis)
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
sf <- read.csv("../data/Robjects/SizeFactors.csv")
pD <- left_join(pD,sf)

#plot(log10(colSums(m))~log10(pD$sf),main="Library Size versus Size Factors",pch=19)
#
#Normalize
m <- log2(t(t(m)/pD$sf)+1)

tsn <- read.csv("../data/Robjects/tSNE.csv")
pD <- left_join(pD,tsn)

ggplot(pD, aes(x=tSNE1, y=tSNE2)) +
    geom_hex(bins=60) +
    facet_wrap(Tissue~Replicate) +
    scale_fill_viridis()

# ggplot(pD, aes(x=tSNE1, y=tSNE2,color=Condition)) +
#     geom_point() +
#     facet_wrap(~Replicate)

# Highly variable genes
hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
hvg <- hvg[rownames(hvg) %in% fD$id[fD$KeepForHvg],]
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

#clustering 
library(igraph)
igr <- buildSNNGraph(m[hvg,],k=15)
clus <- cluster_louvain(igr)
cs <- clus$membership
pD$Cluster <- paste0("C",cs)

ggplot(pD,aes(x=tSNE1, y=tSNE2, color=Cluster)) +
    geom_point()
