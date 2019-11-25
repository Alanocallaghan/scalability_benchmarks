#!/usr/bin/env Rscript

library("here")
source(here("packrat/init.R"))
library("BASiCS")
library("Scalability")

args <- commandArgs(trailingOnly = TRUE)
print(args)
source(here("scripts/chain-scripts/benchmark_code.R"))

set.seed(42)

data <- readRDS(here("data", paste0(args[[1]], ".rds")))
dir <- args[[3]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

frac <- as.numeric(args[[2]])
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
  WithSpikes = TRUE,
  N = 20000,
  Thin = 10,
  Burn = 10000
)
config <- list(
  data = "zeisel",
  chains = 1,
  by = "gene",
  seed = 42,
  proportion_retained = as.numeric(args[[2]])
)
saveRDS(chain, file = file.path(dir, "chains.rds"))
saveRDS(NULL, file = file.path(dir, "time.rds"))
saveRDS(config, file = file.path(dir, "config.rds"))
