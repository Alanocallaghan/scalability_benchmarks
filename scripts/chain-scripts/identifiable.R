args <- commandArgs(trailingOnly = TRUE)
print(args)


library("packrat")

D <- readRDS(paste0("data/", args[[1]], ".rds"))

N <- 20000
Thin <- 10
Burn <- 10000

N <- 10
Thin <- 2
Burn <- 4


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
dir <- paste0(args[[2]], "/", args[[1]])
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(id, file.path(dir, "id.rds"))
