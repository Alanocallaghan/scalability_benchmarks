options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")
library("coda")

theme_set(theme_bw())
source(here("scripts/analysis-scripts/functions.R"))


source(here("scripts/analysis-scripts/identifiability.R"))
source(here("scripts/analysis-scripts/batchinfo.R"))


advi_files <- list.files("outputs/advi", full.names = TRUE)
advi_triplets <- file2triplets(advi_files)
advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
advi_triplets <- lapply(advi_triplets, function(x) x[-3])
advi_df <- read_triplets(advi_triplets)


source(here("scripts/analysis-scripts/elbo_plots.R"))


datasets <- unique(advi_df[["data"]])
data_dims <- vapply(
  datasets,
  function(x) {
    suppressMessages(
      dim(
        readRDS(paste0("data/", x, ".rds"))
      )
    )
  },
  FUN.VALUE = numeric(2)
)
data_dims <- as.data.frame(t(data_dims))
colnames(data_dims) <- c("nGenes", "nCells")
data_dims[["data"]] <- datasets

source(here("scripts/analysis-scripts/downsampling.R"))


dc_files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)

file_df <- rbind(advi_df, dc_df)
df <- merge(file_df, data_dims)

source(here("scripts/analysis-scripts/time_plot.R"))
references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

source(here("scripts/analysis-scripts/de_on_table.R"))
source(here("scripts/analysis-scripts/chain_plots.R"))


source(here("scripts/analysis-scripts/ess.R"))
source(here("scripts/analysis-scripts/hpd.R"))
