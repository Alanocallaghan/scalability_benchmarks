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

sce <- readRDS("rdata/chen.rds")

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


dir.create("outputs/fix_nu", showWarnings = FALSE)
saveRDS(fit_fix, file.path(args[["output"]], "chen-fix.rds"))
saveRDS(fit_var, file.path(args[["output"]], "chen-var.rds"))

droplet_sce <- readRDS("rdata/ibarra-soria.rds")
ind_presom <- colData(droplet_sce)[["Cell_type"]] == "PSM"
ind_som <- colData(droplet_sce)[["Cell_type"]] == "SM"


PSM_Data <- droplet_sce[, ind_presom]
SM_Data <- droplet_sce[, ind_som]

fit_fix <- BASiCS_MCMC(
    PSM_Data,
    PrintProgress = FALSE,
    FixNu = TRUE,
    WithSpikes = FALSE,
    Regression = TRUE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)

fit_var <- BASiCS_MCMC(
    PSM_Data,
    PrintProgress = FALSE,
    FixNu = FALSE,
    WithSpikes = FALSE,
    Regression = TRUE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)
saveRDS(fit_fix, file.path(args[["output"]], "ibarra-presom-fix.rds"))
saveRDS(fit_var, file.path(args[["output"]], "ibarra-presom-var.rds"))

fit_fix <- BASiCS_MCMC(
    SM_Data,
    PrintProgress = FALSE,
    FixNu = TRUE,
    WithSpikes = FALSE,
    Regression = TRUE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)

fit_var <- BASiCS_MCMC(
    SM_Data,
    PrintProgress = FALSE,
    FixNu = FALSE,
    WithSpikes = FALSE,
    Regression = TRUE,
    N = args[["iterations"]],
    Thin = max((args[["iterations"]] / 2) / 1000, 2),
    Burn = max(args[["iterations"]] / 2, 4)
)
saveRDS(fit_fix, file.path(args[["output"]], "ibarra-som-fix.rds"))
saveRDS(fit_var, file.path(args[["output"]], "ibarra-som-var.rds"))
