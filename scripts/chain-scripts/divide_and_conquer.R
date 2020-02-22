#!/usr/bin/env Rscript

library("here")
devtools::load_all("../BASiCS")
# library("Scalability")
# devtools::load_all("../Scalability")
library("future")
plan("multicore")
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)

source(here("scripts/chain-scripts/benchmark_code.R"))


data <- readRDS(here("data", paste0(args[[1]], ".rds")))
dir <- args[[5]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)


data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[[1]],
  SubsetBy = args[[4]],
  NSubsets = args[[2]],
  Seed = args[[3]],
  Regression = TRUE,
  MinGenesPerRBF = 100,
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
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(data[["config"]], file = file.path(dir, "config.rds"))


# grid <- read.delim("data/divide_and_conquer_grid.txt", header=FALSE)
# source(here("scripts/chain-scripts/benchmark_code.R"))
# out <- lapply(1:nrow(grid), function(i) {
#   args <- grid[i, ]
#   data <- readRDS(here("data", paste0(args[[1]], ".rds")))
#   divide_and_conquer_benchmark(
#     Data = data,
#     DataName = args[[1]],
#     SubsetBy = args[[4]],
#     NSubsets = args[[2]],
#     Seed = args[[3]],
#     Regression = TRUE,
#     MinGenesPerRBF = 100,
#     N = 20,
#     Thin = 4,
#     Burn = 8
#   )
# })
