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
parser$add_argument("-s", "--seed")
parser$add_argument("-f", "--fraction")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


source(here("src/chain-scripts/benchmark_code.R"))


data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

set.seed(as.numeric(args[["seed"]]))
frac <- as.numeric(args[["fraction"]])
data <- data[, sample(ncol(counts(data)), floor(ncol(counts(data)) * frac))]
if (length(altExpNames(data))) {
  spikes <- altExp(data, "spike-ins")
  ind_keep_spike <- rowSums(assay(spikes)) != 0
  metadata(data)$SpikeInput <- metadata(data)$SpikeInput[ind_keep_spike, ]
  altExp(data, "spike-ins") <- spikes[ind_keep_spike, ]
}

data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[["data"]],
  SubsetBy = "gene",
  NSubsets = 16,
  Seed = args[["seed"]],
  Regression = TRUE,
  # Verbose = FALSE,
  N = 8,
  Thin = 2,
  Burn = 4
  # N = 20000,
  # Thin = 10,
  # Burn = 10000
)
cfg <- data[["config"]]
cfg$proportion_retained <- as.numeric(args[["fraction"]])
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(cfg, file = file.path(dir, "config.rds"))
