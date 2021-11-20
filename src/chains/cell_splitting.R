#!/usr/bin/env Rscript

## TODO: UNUSED

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
  library("SingleCellExperiment")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


set.seed(args[["seed"]])
sce <- readRDS(sprintf("rdata/%s.rds", args[["data"]]))
fit <- BASiCS_MCMC(
  sce,
  SubsetBy = "cell",
  NSubsets = 4,
  Regression = TRUE,
  PrintProgress = FALSE,
  WithSpikes = "spike-ins" %in% altExpNames(sce),
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4)
)
dir.create("outputs/cell_splitting/", showWarnings = FALSE)
saveRDS(fit,
  sprintf("outputs/cell_splitting/%s_%d.rds", args[["data"]], args[["seed"]])
)
