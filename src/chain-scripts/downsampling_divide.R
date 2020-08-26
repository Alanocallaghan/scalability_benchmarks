#!/usr/bin/env Rscript
if (!require("argparse")) {
    install.packages("argparse")
}
suppressPackageStartupMessages({
  library("argparse")
  library("here")
  # library("BASiCS")
  devtools::load_all("../BASiCS")
  library("future")
  devtools::load_all("../Scalability")
  plan("multicore")
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

counts <- counts(data)
counts[] <- apply(
  counts,
  2,
  function(col) {
    rbinom(length(col), col, as.numeric(args[["fraction"]]))
  }
)
counts(data) <- counts


data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[["data"]],
  SubsetBy = "gene",
  NSubsets = 16,
  Seed = args[["seed"]],
  Regression = TRUE,
  N = 8,
  Thin = 2,
  Burn = 4
  # N = 20000,
  # Thin = 10,
  # Burn = 10000
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
cfg$downsample_rate <- as.numeric(args[["fraction"]])
saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(cfg, file = file.path(dir, "config.rds"))
