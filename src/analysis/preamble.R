suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggbeeswarm")
  library("here")
  library("BASiCS")
  library("coda")
  library("viridis")
})

source(here("src/analysis/functions.R"))
theme_set(theme_bw())

datasets <- c("buettner", "chen", "tung", "zeisel")
data_dims <- vapply(
    datasets,
    function(x) {
        suppressMessages(
            dim(
                readRDS(paste0("rdata/", x, ".rds"))
            )
        )
    },
    FUN.VALUE = numeric(2)
)
data_dims <- as.data.frame(t(data_dims))
colnames(data_dims) <- c("nGenes", "nCells")
data_dims[["data"]] <- datasets
data_dims[["libsize"]] <- vapply(datasets,
    function(x) {
        median(colSums(counts(readRDS(paste0("rdata/", x, ".rds")))))
    },
    numeric(1)
)

advi_files <- list.files("outputs/advi", full.names = TRUE)
advi_triplets <- file2triplets(advi_files)
advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
advi_triplets <- lapply(advi_triplets, function(x) x[-3])
advi_df <- read_triplets(advi_triplets)

dc_files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)

file_df <- rbind(advi_df, dc_df)
# file_df <- dc_df
df <- merge(file_df, data_dims)

references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)
