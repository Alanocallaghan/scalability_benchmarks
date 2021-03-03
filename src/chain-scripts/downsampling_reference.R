#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-f", "--fraction", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()



data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)


set.seed(42)
counts <- counts(data)
counts[] <- apply(
  counts,
  2,
  function(col) {
    rbinom(length(col), col, as.numeric(args[["fraction"]]))
  }
)
counts(data) <- counts

time <- system.time(
  chain <- BASiCS_MCMC(
    data,
    WithSpikes = length(SingleCellExperiment::altExpNames(data)) > 0,
    Regression = TRUE,
    PrintProgress = FALSE,
    # N = 20000,
    # Thin = 10,
    # Burn = 10000
    N = 8,
    Thin = 2,
    Burn = 4
  )
)
config <- list(
  data = args[["data"]],
  chains = 1,
  by = "gene",
  seed = 42,
  downsample_rate = as.numeric(args[["fraction"]])
)
saveRDS(chain, file = file.path(dir, "chains.rds"))
saveRDS(time, file = file.path(dir, "time.rds"))
saveRDS(config, file = file.path(dir, "config.rds"))
