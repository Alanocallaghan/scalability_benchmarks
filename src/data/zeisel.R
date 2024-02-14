library("scRNAseq")
library("scater")

Data <- scRNAseq::ZeiselBrainData()
Data <- Data[, Data$level1class == "pyramidal CA1"]

libsize_drop <- isOutlier(
    colSums(counts(Data)),
    nmads = 3,
    type = "lower",
    log = TRUE
)
feature_drop <- isOutlier(
    colSums(counts(Data) != 0),
    nmads = 3,
    type = "lower",
    log = TRUE
)
Data <- Data[, !(libsize_drop | feature_drop)]

bio <- counts(Data)[!grepl("ERCC", rownames(Data)), , drop = FALSE]
ind_expressed <- rowMeans(bio) >= 1 & rowMeans(bio != 0) > 0.5

ind_expressed_ercc <- c(
    which(ind_expressed),
    grep("ERCC", rownames(Data))
)
Data <- Data[ind_expressed_ercc, ]

Data$BatchInfo <- Data$group
altExp(Data, "repeat") <- NULL
saveRDS(Data, file = "rdata/zeisel.rds")
