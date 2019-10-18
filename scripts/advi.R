#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

library("packrat")
source("packrat/init.R")
library("here")
library("BASiCS")
library("Scalability")


set.seed(args[[2]])
data <- readRDS(here("data", paste0(args[[1]], ".rds")))
dir <- args[[3]]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

with_spikes <- as.logical(length(altExpNames(data)))
time <- system.time(
  chain <- Scalability:::BASiCS_stan(
    data,
    WithSpikes = with_spikes,
    Regression = TRUE
#    ,
#    iter = 10,
#    tol_rel_obj = 1,
#    eta = 0.1,
#    eval_elbo = 2,
#    adapt_engaged = FALSE
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
  data = args[[1]],
  seed = NA
)

saveRDS(time, file.path(dir, "time.rds"))
saveRDS(config, file.path(dir, "config.rds"))
saveRDS(chain, file.path(dir, "chain.rds"))
