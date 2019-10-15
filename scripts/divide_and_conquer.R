#!/usr/bin/env Rscript

library("packrat")
library("here")
library("BASiCS")
library("Scalability")

args <- commandArgs(trailingOnly = TRUE)
source(here("scripts/benchmark_code.R"))


data <- get(args[[1]])
dir <- args[[5]]
dir.create(dir, recursive = TRUE)
print(dir)

data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[[1]],
  SubsetBy = "gene",
  NSubsets = args[[2]],
  Seed = args[[3]],
  Regression = TRUE,
  Verbose = FALSE,
  N = 10,
  Thin = 2,
  Burn = 4
)
chains <- data[["chain"]]
config <- data[["config"]]
if (length(chains) == 1) {
  collapsed <- chains
} else {
  collapsed <- combine_subposteriors(
    chains,
    subset_by = config[["by"]],
    method = "pie",
    weight_method = "n_weight",
    mc.cores = 1
  )
}
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(data[["config"]], file = file.path(dir, "config.rds"))
