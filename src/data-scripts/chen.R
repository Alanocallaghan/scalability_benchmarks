library("scater")
library("scran")


sce <- readRDS("downloads/chen.rds")

sce <- sce[, sce$cell_type1 == "Astro"]

# hist(log10(colSums(counts(sce))), breaks = "FD")
# hist(log10(rowMeans(counts(sce))), breaks = "FD")
# hist(rowMeans(counts(sce) != 0), breaks = "FD")


sce <- sce[rowMeans(counts(sce)) > 0.1, ]
saveRDS(sce, "data/chen.rds")
