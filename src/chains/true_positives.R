#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-n", "--nsubsets", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

# print(dput(args))
# stop()
set.seed(args[["seed"]])

droplet_sce <- readRDS("rdata/ibarra-soria.rds")

# Presomitic mesoderm
ind_presom <- colData(droplet_sce)[["Cell_type"]] == "PSM"
ind_som <- colData(droplet_sce)[["Cell_type"]] == "SM"

PSM_Data <- droplet_sce[, ind_presom]
SM_Data <- droplet_sce[, ind_som]

# Presomitic mesoderm cells
PSM_MCMC <- BASiCS_MCMC(
    PSM_Data,
    Regression = TRUE,
    SubsetBy = "gene",
    NSubsets = args[["nsubsets"]],
    WithSpikes = FALSE,
    PrintProgress = FALSE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)

# Somitic mesoderm cells
SM_MCMC <- BASiCS_MCMC(
    SM_Data,
    Regression = TRUE,
    SubsetBy = "gene",
    NSubsets = args[["nsubsets"]],
    WithSpikes = FALSE,
    PrintProgress = FALSE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)



test <- BASiCS_TestDE(PSM_MCMC, SM_MCMC)

out <- list(
    test = test,
    mcmc = list(sm = SM_MCMC, psm = PSM_MCMC)
)
dir.create("outputs/true-positives/", showWarnings = FALSE)

saveRDS(out, args[["output"]])
