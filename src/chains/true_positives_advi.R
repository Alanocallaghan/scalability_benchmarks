#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCStan")
})
parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

set.seed(args[["seed"]])

droplet_sce <- readRDS("rdata/ibarra-soria.rds")

# Presomitic mesoderm
ind_presom <- colData(droplet_sce)[["Cell_type"]] == "PSM"
ind_som <- colData(droplet_sce)[["Cell_type"]] == "SM"

PSM_Data <- droplet_sce[, ind_presom]
SM_Data <- droplet_sce[, ind_som]

# Presomitic mesoderm cells
PSM_MCMC <- BASiCStan(
    PSM_Data,
    Regression = TRUE,
    WithSpikes = FALSE
)

# Somitic mesoderm cells
SM_MCMC <- BASiCStan(
    SM_Data,
    Regression = TRUE,
    WithSpikes = FALSE
)

test <- BASiCS_TestDE(PSM_MCMC, SM_MCMC)

out <- list(
    test = test,
    mcmc = list(sm = SM_MCMC, psm = PSM_MCMC)
)

dir.create("outputs/true-positives/", showWarnings = FALSE)
saveRDS(out, args[["output"]])
