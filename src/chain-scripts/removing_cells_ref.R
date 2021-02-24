#!/usr/bin/env Rscript
if (!require("argparse")) {
    install.packages("argparse")
}
suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-f", "--fraction")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

source(here("src/chain-scripts/benchmark_code.R"))

set.seed(42)
data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

frac <- as.numeric(args[["fraction"]])
data <- data[, sample(ncol(counts(data)), floor(ncol(counts(data)) * frac))]
if (length(altExpNames(data))) {
  spikes <- altExp(data, "spike-ins")
  ind_keep_spike <- rowSums(assay(spikes)) != 0
  metadata(data)$SpikeInput <- metadata(data)$SpikeInput[ind_keep_spike, ]
  altExp(data, "spike-ins") <- spikes[ind_keep_spike, ]
}


chain <- BASiCS_MCMC(
  data,
  Regression = TRUE,
  WithSpikes = as.logical(length(altExpNames(data))),
  N = 8,
  Thin = 2,
  Burn = 4
  # N = 20000,
  # Thin = 10,
  # Burn = 10000
)
config <- list(
  data = args[["data"]],
  chains = 1,
  by = "gene",
  seed = 42,
  proportion_retained = as.numeric(args[["fraction"]])
)
saveRDS(chain, file = file.path(dir, "chains.rds"))
saveRDS(NULL, file = file.path(dir, "time.rds"))
saveRDS(config, file = file.path(dir, "config.rds"))
