#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("here")
    library("BASiCS")
    library("BASiCStan")
})

sce <- readRDS("rdata/tung.rds")

set.seed(42)
time_hmc <- system.time(
    fit_hmc <- BASiCStan(
        sce,
        Method = "sampling",
        chains = 1,
        Verbose = TRUE,
        WithSpikes = TRUE,
        Regression = TRUE
    )
)

time_amwg <- system.time(
    fit_amwg <- BASiCS_MCMC(
        sce,
        PrintProgress = FALSE,
        WithSpikes = TRUE,
        Regression = TRUE,
        N = 40000,
        Thin = 10,
        Burn = 20000
    )
)

out <- list(
    chains = list(
        hmc = fit_hmc,
        amwg = fit_amwg
    ),
    time = list(
        hmc = time_hmc,
        amwg = time_amwg
    )
)
saveRDS(out, "outputs/hmc_vs_amwg.rds")
