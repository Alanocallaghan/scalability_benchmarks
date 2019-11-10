#!/usr/bin/env Rscript

library("here")
source(here("packrat/init.R"))
library("BASiCS")
library("Scalability")

args <- commandArgs(trailingOnly = TRUE)

data <- readRDS(here("data", paste0(args[[1]], ".rds")))
dir <- args[[3]]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)


counts <- counts(data)
counts[] <- apply(
  counts,
  2,
  function(col) {
    rbinom(length(col), col, as.numeric(args[[2]]))
  }
)
counts(data) <- counts
set.seed(42)

time <- system.time(
  chain <- BASiCS_MCMC(
    data,
    WithSpikes = TRUE,
    Regression = TRUE,
    N = 20000,
    Thin = 10,
    Burn = 10000
  )
)
config <- list(
  data = "zeisel",
  chains = 1,
  by = "gene",
  seed = 42,
  downsample_rate = as.numeric(args[[2]])
)
saveRDS(chain, file = file.path(dir, "chains.rds"))
saveRDS(time, file = file.path(dir, "time.rds"))
saveRDS(config, file = file.path(dir, "config.rds"))
