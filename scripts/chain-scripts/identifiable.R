args <- commandArgs(trailingOnly = TRUE)
print(args)


library("packrat")

data <- readRDS(paste0("data/", args[[1]], ".rds"))

N <- 20000
Thin <- 10
Burn <- 10000

id <- with_extlib(
  c(
    "SingleCellExperiment",
    "BASiCS"
  ), 
    BASiCS_MCMC(
      data,
      N = N,
      Thin = Thin,
      Burn = Burn,
      WithSpikes = FALSE,
      Regression = TRUE
    )
  }
)

saveRDS(id, paste0(args[[2]], "/", args[[1]], "/id.rds"))
