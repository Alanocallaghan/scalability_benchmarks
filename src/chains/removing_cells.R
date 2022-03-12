#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-f", "--fraction", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

source(here("src/chains/benchmark_code.R"))


data <- readRDS(here("rdata", paste0(args[["data"]], ".rds")))
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
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)
cfg <- data[["config"]]
cfg$proportion_retained <- as.numeric(args[["fraction"]])
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(cfg, file = file.path(dir, "config.rds"))
