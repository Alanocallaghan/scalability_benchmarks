library("BASiCS")
args <- commandArgs(trailingOnly = TRUE)
print(args)

non_id <- BASiCS_MCMC(
  data,
  N = N,
  Thin = Thin,
  Burn = Burn,
  WithSpikes = FALSE,
  Regression = TRUE
)
saveRDS(non_id, paste0(args[[2]], "/", args[[1]], "/non_id.rds"))
