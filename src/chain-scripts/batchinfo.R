#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


data <- readRDS(paste0("data/", args[["data"]], ".rds"))

# N <- 20000
# Thin <- 10
# Burn <- 10000
N <- 8
Thin <- 2
Burn <- 4
library("BASiCS")

bi <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  PrintProgress = FALSE,
  WithSpikes = FALSE,
  Regression = TRUE
)

data@colData$BatchInfo <- 1

non_bi <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  PrintProgress = FALSE,
  WithSpikes = FALSE,
  Regression = TRUE
)

dir <- args[["output"]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(bi, file.path(dir, "batch.rds"))
saveRDS(non_bi, file.path(dir, "/nobatch.rds"))
