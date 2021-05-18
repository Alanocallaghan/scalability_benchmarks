#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-n", "--nsubsets", type = "double")
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-b", "--subsetby")
parser$add_argument("-o", "--output")
parser$add_argument("-i", "--iterations", type = "double")
args <- parser$parse_args()

source(here("src/chain-scripts/benchmark_code.R"))

data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)


data <- divide_and_conquer_benchmark(
  Data = data,
  DataName = args[["data"]],
  SubsetBy = args[["subsetby"]],
  NSubsets = args[["nsubsets"]],
  Seed = args[["seed"]],
  Regression = TRUE,
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4)
)
chains <- data[["chain"]]
config <- data[["config"]]

saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(data[["config"]], file = file.path(dir, "config.rds"))
