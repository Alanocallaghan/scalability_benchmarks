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
args <- parser$parse_args()

args <- list(data = "zeisel", nsubsets = 1, output = "outputs/divide_and_conquer/data-zeisel_nsubsets-1_seed-28_by-gene", 
    seed = 28, subsetby = "gene")


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
  N = 20000,
  Thin = 10,
  Burn = 10000
)
chains <- data[["chain"]]
config <- data[["config"]]

saveRDS(data[["chain"]], file = file.path(dir, "chains.rds"))
saveRDS(data[["time"]], file = file.path(dir, "time.rds"))
saveRDS(data[["config"]], file = file.path(dir, "config.rds"))
