library("scran")
library("scater")
library("BASiCS")


if (!file.exists("downloads/rawCounts.tsv")) {
  website <- "https://www.ebi.ac.uk/"
  folder <- "arrayexpress/files/E-MTAB-6153/"
  file <- "E-MTAB-6153.processed.2.zip"
  download.file(
    paste0(website, folder, file),
    destfile = "downloads/rawCounts.zip"
  )
  unzip(zipfile = "downloads/rawCounts.zip", exdir = "downloads")
  file.remove("downloads/rawCounts.zip")
}
rawCounts <- read.delim("downloads/rawCounts.tsv", header = TRUE)


if (!file.exists("downloads/cellAnnotation.tsv")) {
  website <- "https://www.ebi.ac.uk/"
  folder <- "arrayexpress/files/E-MTAB-6153/"
  file <- "E-MTAB-6153.processed.3.zip"
  download.file(
    paste0(website, folder, file),
    destfile = "cluster_labels.zip"
  )
  unzip(zipfile = "cluster_labels.zip", exdir = "downloads")
}

cluster_labels <- read.table("downloads/cellAnnotation.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)


cluster_labels[["Cell_type"]] <- cluster_labels$cellType
cluster_labels[["Cell_type"]] <- sub(
  "presomiticMesoderm",
  "PSM",
  cluster_labels[["Cell_type"]]
)
cluster_labels[["Cell_type"]] <- sub("somiticMesoderm",
  "SM", 
  cluster_labels[["Cell_type"]]
)

ind_som <- which(cluster_labels[["Cell_type"]] == "PSM" |
  cluster_labels[["Cell_type"]] == "SM")
rawCounts <- rawCounts[, ind_som]
cluster_labels <- cluster_labels[ind_som, ]


droplet_sce <- SingleCellExperiment(
  assays = list(counts = as(as.matrix(rawCounts), "dgCMatrix"))
)
rm(rawCounts)
colData(droplet_sce) <- DataFrame(
  cluster_labels,
  subCellType = sub("_.*", "", cluster_labels$cell)
)
ind_expressed <- Matrix::rowMeans(counts(droplet_sce)) > 0.1
droplet_sce <- droplet_sce[ind_expressed, ]

droplet_sce <- computeSumFactors(droplet_sce,
  clusters = colData(droplet_sce)$subCellType
)
droplet_sce <- logNormCounts(droplet_sce)
droplet_sce <- runPCA(droplet_sce)

# Cell types identified by clustering
# plotReducedDim(droplet_sce, dimred = "PCA", colour_by = "subCellType") +
#   scale_fill_manual(name = "cellType", values = c("coral4", "steelblue", "limegreen"))

ind_retain <- reducedDims(droplet_sce)$PCA[, 2] > -5 &
  colData(droplet_sce)$subCellType != "presomiticMesoderm.b"
droplet_sce <- droplet_sce[, ind_retain]

droplet_sce$batch <- as.character(round(droplet_sce$sample))
