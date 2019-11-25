#!/usr/bin/env Rscript

library("here")
source(here("packrat/init.R"))
library("BASiCS")
library("Scalability")

args <- commandArgs(trailingOnly = TRUE)

source(here("scripts/chain-scripts/benchmark_code.R"))


data <- readRDS(here("data", paste0(args[[1]], ".rds")))
dir <- args[[4]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

set.seed(as.numeric(args[[3]]))
frac <- as.numeric(args[[2]])
data <- data[, sample(ncol(counts(data)), floor(ncol(counts(data)) * frac))]
if (length(altExpNames(data))) {
  spikes <- altExp(data, "spike-ins")
  ind_keep_spike <- rowSums(assay(spikes)) != 0
  metadata(data)$SpikeInput <- metadata(data)$SpikeInput[ind_keep_spike, ]
  altExp(data, "spike-ins") <- spikes[ind_keep_spike, ]
}

data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[[1]],
  SubsetBy = "gene",
  NSubsets = 16,
  Seed = args[[3]],
  Regression = TRUE,
  Verbose = FALSE,
  N = 20000,
  Thin = 10,
  Burn = 10000
)
cfg <- data[["config"]]
cfg$proportion_retained <- as.numeric(args[[2]])
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(cfg, file = file.path(dir, "config.rds"))
