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
  WithSpikes = FALSE,
  Regression = TRUE
)

data@colData$BatchInfo <- 1

non_bi <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  WithSpikes = FALSE,
  Regression = TRUE
)

dir <- args[[2]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(bi, file.path(dir, "batch.rds"))
saveRDS(non_bi, file.path(dir, "/nobatch.rds"))
