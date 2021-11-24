library("scater")
library("scRNAseq")
library("BASiCS")

bo <- BuettnerESCData()
b <- bo[, bo$phase == "G2M"]
libsize_drop <- isOutlier(
  Matrix::colSums(counts(b)),
  nmads = 3,
  type = "lower",
  log = TRUE
)
feature_drop <- isOutlier(
  Matrix::colSums(counts(b) != 0),
  nmads = 3,
  type = "lower",
  log = TRUE
)
ind_drop <- !(libsize_drop | feature_drop)
b <- b[, ind_drop]
ind_expressed <- (Matrix::rowMeans(counts(b)) >= 100 & 
  Matrix::rowMeans(counts(b) != 0) > 0.5)
b <- b[ind_expressed, ]

erc <- altExp(b, "ERCC")
altExp(b, "ERCC") <- NULL
nospike <- rowSums(counts(erc)) == 0
erc <- erc[!nospike, ]
rowData(erc) <- data.frame(id = rownames(erc), rowData(erc)[, "molecules"])
altExp(b, "spike-ins") <- erc

saveRDS(sce, "rdata/buettner.rds")
