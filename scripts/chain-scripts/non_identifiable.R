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
dir <- paste0(args[[2]], "/", args[[1]])
dir.create(dir, showWarnings = FALSE)
saveRDS(non_id, file.path(dir, "non_id.rds"))
