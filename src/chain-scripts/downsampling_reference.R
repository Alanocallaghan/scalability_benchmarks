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
parser$add_argument("-f", "--fraction")
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
    N = 8,
    Thin = 2,
    Burn = 4
    # N = 20000,
    # Thin = 10,
    # Burn = 10000
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
