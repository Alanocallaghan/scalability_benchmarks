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
  elbo <- capture.output(
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


## Convergence first
parse_elbo <- function(c) {
  c <- gsub("Chain 1:\\s+", "", c)
  elbo <- c[(grep("Begin stochastic", c) + 1):(grep("Drawing", c) - 2)]
  elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
    elbo[-c(1, grep("CONVERGED", elbo))],
    "NOTCONVERGED")
  elbo <- gsub("(MEDIAN )?ELBO CONVERGED", "CONVERGED", elbo)
  elbo <- strsplit(elbo, "\\s+")
  elbo <- do.call(rbind, elbo)
  colnames(elbo) <- elbo[1, ]
  elbo <- elbo[-1, ]
  elbo <- as.data.frame(elbo, stringsAsFactors=FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}


saveRDS(parse_elbo(elbo), file.path(dir, "elbo.rds"))
saveRDS(time, file.path(dir, "time.rds"))
saveRDS(config, file.path(dir, "config.rds"))
saveRDS(chain, file.path(dir, "chain.rds"))
