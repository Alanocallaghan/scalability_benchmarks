#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("argparse")
    library("here")
    library("ggplot2")
    library("BASiCS")
    library("viridis")
    library("ggpointdensity")
    library("cowplot")
})
parser <- ArgumentParser()
parser$add_argument("-f", "--fix")
parser$add_argument("-v", "--var")
parser$add_argument("-d", "--dataset")
args <- parser$parse_args()

# if (is.null(args[["input"]])) {
#     args[["input"]] <- "outputs/fix_nu"
# }
source("src/analysis/functions.R")

theme_set(theme_bw())

fit_fix <- readRDS(args[["fix"]])
fit_var <- readRDS(args[["var"]])

summary_var <- Summary(fit_var)
summary_fix <- Summary(fit_fix)

ord <- order(summary_fix@parameters$mu[, "median"])


g1 <- plot_hpds(summary_var, summary_fix, "mu", ord)
g2 <- plot_hpds(summary_var, summary_fix, "delta", ord)
g3 <- plot_hpds(summary_var, summary_fix, "epsilon", ord)

combined <- plot_with_legend_below(g1, g2, g3)
ggsave(sprintf("figs/fixnu-%s.pdf", args[["dataset"]]), width = 8, height = 3)


d1 <- plot_hpd_diff(summary_var, summary_fix, "mu", ord)
d2 <- plot_hpd_diff(summary_var, summary_fix, "delta", ord)
d3 <- plot_hpd_diff(summary_var, summary_fix, "epsilon", ord)

dop <- plot_grid(d1, d2, d3, labels = "AUTO", nrow = 1)
ggsave(sprintf("figs/fixnu-diff-%s.pdf", args[["dataset"]]), width = 8, height = 3)
