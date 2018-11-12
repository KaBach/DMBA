
## This script extracts a cleaned count matrix from the 10X output
## The script is largely based on:
## https://github.com/MarioniLab/DropletUtils

# First remove reads derived from index swapping
library(DropletUtils)
samples <-  c("../../CellRangerAnalysis/SampleA10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleB10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleC10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleD10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleE10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleF10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleG10/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleH10/outs/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
		     barcode.length=bc.length,
		     get.swapped=get.swapped,
		     min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("A10","B10","C10","D10","E10","F10","G10","H10")

finalCounts <- list()
# bpparam <- MulticoreParam(workers=3)

for (i in seq_along(cleanedCounts)) {
    m <- cleanedCounts[[i]]

    # test for empty droplets
    out <- emptyDrops(m)#,BPPARAM=bpparam)

    # set FDR threshold
    nonempty <- out$FDR < 0.01

    #remove NAs from NA FDR
    nonempty[is.na(nonempty)] <- FALSE

    # cant use logic vector to subset dgCMatrix?
    nonempty <- rownames(out[nonempty,])

    # Clean matrix
    m.clean <- m[,nonempty]
    finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../data/Robjects/CountMatrix.rds")
