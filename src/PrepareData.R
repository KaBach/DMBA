# Script to prepare cellranger data for downstream analysis
library(plyr)
library(dplyr)
library(reshape2)

# ---- ReadData ----

# Read in output dropletUtils package
cDat <- readRDS("../data/Robjects/CountMatrix.rds")
pDat <- data.frame(barcode=colnames(cDat))
fDat <- read.table("../data/CellRangerData/SampleA10/outs/filtered_gene_bc_matrices/mm10_tg/genes.tsv",stringsAsFactors=FALSE)
colnames(fDat) <- c("id","symbol")
rownames(fDat) <- fDat$id
fDat <- fDat[rownames(cDat),]

#reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
fDat <- fDat[keep,]

# ---- Formatting ----

# conversion into sample names
conv <- read.csv("../data/misc/SLX-16967.H2T5LBBXY.s_6.contents.csv", stringsAsFactors=FALSE)
pDat$barcode <- as.character(pDat$barcode)
conv$Sample.name <- gsub("-","",conv$Sample.name) 
conv$Barcode <- gsub("SIGA","",conv$Barcode)

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,nchar(barcode))) %>%
        mutate(SampleID=mapvalues(SeqID,conv$Barcode,
				  conv$Sample.name)) %>%
	mutate(Tissue=mapvalues(SeqID,conv$Barcode,
				conv$Tissue)) %>%
	mutate(Red=mapvalues(SeqID,conv$Barcode,
			     conv$Red)) %>%
	mutate(Condition=paste(Tissue,Red,sep=".")) %>%
	mutate(Replicate=as.factor(gsub("[A-Z]","",SampleID)))

# Add more info to the feature Data
mitoGenes <- read.table("../data/misc/MitoGenes.txt")
tfCheck <- read.table("../data/misc/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")

# Ribosmal genes
library(biomaRt)
gos <- c("GO:0003735","G0:0005840","GO:0015935","GO:0015934") # constituents of ribosome
ensembl  <- useMart("ensembl",dataset="mmusculus_gene_ensembl") 
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
gene.data <- gene.data[gene.data$go_id %in% gos,]

fDat$Ribosomal <- fDat$id %in% gene.data$ensembl_gene_id

# Mitochondrial
fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1

# TranscriptionFactor
fDat$TranscriptionFactor <- fDat$id %in% tfCheck$ensembl_gene_id

# Keep for hvg

fDat$KeepForHvg <- !(fDat$Ribosomal | fDat$Mitochondrial)


# Save data
stopifnot(identical(rownames(fDat),as.character(rownames(cDat))) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/ExpressionList.rds")
