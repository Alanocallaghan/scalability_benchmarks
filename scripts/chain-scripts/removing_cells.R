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

frac <- as.numeric(args[[2]])

sample(ncol(counts), floor(ncol(counts) * frac))

counts <- counts(data)

ref <- BASiCS_MCMC(
  data,
  N = 20000,
  Thin = 10,
  Burn = 10000
)

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
chains <- data[["chain"]]
config <- data[["config"]]
if (length(chains) == 1) {
  collapsed <- chains
} else {
  collapsed <- Scalability:::combine_subposteriors(
    chains,
    subset_by = config[["by"]],
    method = "pie",
    weight_method = "n_weight",
    mc.cores = 1
  )
}
cfg <- data[["config"]]
cfg$downsample_rate <- as.numeric(args[[2]])
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(cfg, file = file.path(dir, "config.rds"))
