args <- commandArgs(trailingOnly = TRUE)

library("packrat")

data <- readRDS(paste0("data/", args[[1]], ".rds"))

N <- 20000
Thin <- 10
Burn <- 10000
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

dir <- paste0(args[[2]], "/", args[[1]], "/")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(bi, file.path(dir, "batch.rds"))
saveRDS(non_bi, file.path(dir, "/nobatch.rds"))
