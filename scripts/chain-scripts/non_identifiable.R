library("BASiCS")
args <- commandArgs(trailingOnly = TRUE)
print(args)

D <- readRDS(paste0("data/", args[[1]], ".rds"))


non_id <- BASiCS_MCMC(
  D,
  N = 20000,
  Thin = 10,
  Burn = 10000,
  WithSpikes = FALSE,
  Regression = TRUE
)
dir <- paste0(args[[2]], "/", args[[1]])
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(non_id, file.path(dir, "non_id.rds"))
