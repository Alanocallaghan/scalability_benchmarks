library("edgeR")
library("scran")
library("scater")
library("Seurat")


# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3

options(stringAsFactors = FALSE)

##https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-018-1578-4
## This tells us ENSG00000105374, a marker for nk cells
system("curl https://static-content.springer.com/esm/art%3A10.1186%2Fs12967-018-1578-4/MediaObjects/12967_2018_1578_MOESM2_ESM.xlsx > data/pbmc_markers.xlsx")


markers <- readxl::read_xlsx("data/pbmc_markers.xlsx", skip = 2)

x10 <- read10X(
  "downloads/5k_pbmc_protein_v3_filtered_feature_bc_matrix/matrix.mtx.gz",
  "downloads/5k_pbmc_protein_v3_filtered_feature_bc_matrix/features.tsv.gz",
  "downloads/5k_pbmc_protein_v3_filtered_feature_bc_matrix/barcodes.tsv.gz"
)

sce <- SingleCellExperiment(
  assays = list(counts = x10$counts),
  colData = x10$samples,
  rowData = x10$genes
)


clusters <- quickCluster(sce, min.size = 100)
sce <- computeSumFactors(sce, cluster = clusters)
sce <- normalize(sce)



markers <- markers[markers$Ensembl %in% rownames(sce), ]

# Heatmap(logcounts(d)[markers$Ensembl, ], row_split = markers$Subpopulation)

## Subset based on ENSG00000105374 (as above)
sce <- sce[, logcounts(sce)["ENSG00000105374", ] > 3]


ind_gene <- rowMeans(counts(sce) != 0) > 0.2
ind_cell <- colSums(counts(sce)) > 1000

sce <- sce[ind_gene, ind_cell]
counts(sce) <- as.matrix(counts(sce))
libsize_drop <- isOutlier(
  Matrix::colSums(counts(sce)),
  nmads = 3,
  type = "lower",
  log = TRUE
)
feature_drop <- isOutlier(
  Matrix::colSums(counts(sce) != 0),
  nmads = 3,
  type = "lower",
  log = TRUE
)
ind_expressed <- Matrix::rowMeans(counts(sce)) > 0.1 &
  Matrix::rowMeans(counts(sce) != 0) >= 0.2

ind_drop <- libsize_drop | feature_drop
sce <- sce[ind_expressed, !ind_drop]
colnames(sce) <- paste0("Cell", seq_len(ncol(sce)))

pbmc <- sce
saveRDS(pbmc, "data/pbmc.rds")
# usethis::use_data(pbmc, overwrite = TRUE)
# saveRDS(sce, "datasets/pbmc_nk_data.rds")

