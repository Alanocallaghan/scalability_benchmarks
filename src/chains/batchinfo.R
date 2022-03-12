#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


data <- readRDS(paste0("rdata/", args[["data"]], ".rds"))


bi <- BASiCS_MCMC(
    data,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4),
    PrintProgress = FALSE,
    WithSpikes = FALSE,
    Regression = TRUE
)

data@colData$BatchInfo <- 1

non_bi <- BASiCS_MCMC(
    data,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4),
    PrintProgress = FALSE,
    WithSpikes = FALSE,
    Regression = TRUE
)

dir <- args[["output"]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(bi, file.path(dir, "batch.rds"))
saveRDS(non_bi, file.path(dir, "/nobatch.rds"))
