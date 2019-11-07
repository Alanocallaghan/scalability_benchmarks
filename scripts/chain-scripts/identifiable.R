args <- commandArgs(trailingOnly = TRUE)

library("packrat")

data <- readRDS(paste0("data/", args[[1]], ".rds"))
data@colData$BatchInfo <- sample(2, ncol(data@assays@data$counts), replace = TRUE)

N <- 20000
Thin <- 10
Burn <- 10000

id <- with_extlib("BASiCS", 
  BASiCS_MCMC(
    data,
    N = N,
    Thin = Thin,
    Burn = Burn,
    WithSpikes = FALSE,
    Regression = TRUE
  )
)

library("BASiCS")

non_id <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  WithSpikes = FALSE,
  Regression = TRUE
)

saveRDS(id, paste0(args[[2]], "/", args[[1]], "/id.rds"))
saveRDS(non_id, paste0(args[[2]], "/", args[[1]], "/non_id.rds"))
