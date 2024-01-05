library("scater")
library("scran")

if (!file.exists("downloads/GSE87544_Merged_17samples_14437cells_count.txt")) {
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87544/suppl/GSE87544%5FMerged%5F17samples%5F14437cells%5Fcount%2Etxt%2Egz", 
        destfile = "downloads/GSE87544_Merged_17samples_14437cells_count.txt.gz")
    system("gunzip downloads/GSE87544_Merged_17samples_14437cells_count.txt.gz")
}
if (!file.exists("downloads/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv")) {
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87544/suppl/GSE87544%5F1443737Cells%2ESVM%2Ecluster%2Eidentity%2Erenamed%2Ecsv%2Egz", 
        destfile = "downloads/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz")
    system("gunzip downloads/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz")
}

### https://github.com/hemberg-lab/scRNA.seq.datasets/blob/5cb8c3896f1d70a7143070164b0e4102e225c323/R/chen.R
### DATA
d <- read.table("downloads/GSE87544_Merged_17samples_14437cells_count.txt", header = TRUE)
rownames(d) <- d[,1]
d <- d[,2:ncol(d)]
d <- d[,order(colnames(d))]

### ANNOTATIONS
ann <- read.csv("downloads/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv")
rownames(ann) <- ann[,2]
colnames(ann)[3] <- "cell_type1"
ann <- ann[,3,drop = FALSE]
ann <- ann[order(rownames(ann)), , drop = FALSE]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(d)), colData = ann)
sce <- sce[, sce$cell_type1 == "Astro"]

# hist(log10(colSums(counts(sce))), breaks = "FD")
# hist(log10(rowMeans(counts(sce))), breaks = "FD")
# hist(rowMeans(counts(sce) != 0), breaks = "FD")


sce <- sce[rowMeans(counts(sce)) > 0.1, ]
saveRDS(sce, "rdata/chen.rds")
