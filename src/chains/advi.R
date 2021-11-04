#!/usr/bin/env Rscript

library("argparse")
library("here")
library("BASiCS")
library("BASiCStan")

parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


set.seed(args[["seed"]])
data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

with_spikes <- "spike-ins" %in% altExpNames(data)
time <- system.time(
  elbo <- capture.output(
    chain <- BASiCStan(
      data,
      WithSpikes = with_spikes
    )
  )
)
config <- list(
  chains = NA,
  by = "advi",
  data = args[["data"]],
  seed = args[["seed"]]
)

saveRDS(elbo, file.path(dir, "elbo.rds"))
saveRDS(time, file.path(dir, "time.rds"))
saveRDS(config, file.path(dir, "config.rds"))
saveRDS(chain, file.path(dir, "chain.rds"))
