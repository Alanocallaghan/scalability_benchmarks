args <- commandArgs(trailingOnly = TRUE)

D <- readRDS(paste0("data/", args[[1]], ".rds"))

N <- 20000
Thin <- 10
Burn <- 10000

id <- with_extlib(
  c(
    "SingleCellExperiment",
    "BASiCS"
  ), 
  BASiCS_MCMC(
    D,
    N = N,
    Thin = Thin,
    Burn = Burn,
    WithSpikes = FALSE,
    Regression = TRUE
  )
)
dir <- args[[2]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(id, file.path(dir, "id.rds"))
