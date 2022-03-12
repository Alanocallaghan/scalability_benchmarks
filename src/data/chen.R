library("scater")
library("scran")

if (!file.exists("downloads/chen.rds")) {
    system("curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/chen.rds > downloads/chen.rds")
}

sce <- readRDS("downloads/chen.rds")

sce <- sce[, sce$cell_type1 == "Astro"]

# hist(log10(colSums(counts(sce))), breaks = "FD")
# hist(log10(rowMeans(counts(sce))), breaks = "FD")
# hist(rowMeans(counts(sce) != 0), breaks = "FD")


sce <- sce[rowMeans(counts(sce)) > 0.1, ]
saveRDS(sce, "rdata/chen.rds")
