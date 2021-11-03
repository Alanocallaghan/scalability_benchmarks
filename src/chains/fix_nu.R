#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

sce <- readRDS("data/chen.rds")

fit_fix <- BASiCS_MCMC(
  sce,
  PrintProgress = FALSE,
  FixNu = TRUE,
  WithSpikes = FALSE,
  Regression = TRUE,
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4)
)

fit_var <- BASiCS_MCMC(
  sce,
  PrintProgress = FALSE,
  FixNu = FALSE,
  WithSpikes = FALSE,
  Regression = TRUE,
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4)
)

dir.create(args[["output"]])
saveRDS(fit_fix, file.path(args[["output"]], "fix.rds"))
saveRDS(fit_var, file.path(args[["output"]], "var.rds"))
