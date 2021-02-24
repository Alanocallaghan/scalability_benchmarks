#!/usr/bin/env Rscript
if (!require("argparse")) {
    install.packages("argparse")
}
suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-s", "--seed")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


set.seed(args[["seed"]])
data <- readRDS(here("data", paste0(args[["data"]], ".rds")))
dir <- args[["output"]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

with_spikes <- as.logical(length(altExpNames(data)))
time <- system.time(
  elbo <- capture.output(
    chain <- Scalability:::BASiCS_stan(
      data,
      WithSpikes = with_spikes,
      Regression = TRUE
    )
  )
)
chain <- Scalability:::stan2basics(
  chain, 
  gene_names = rownames(counts(data)),
  cell_names = colnames(data)
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
