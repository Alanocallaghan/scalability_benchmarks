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
parser$add_argument("-i", "--input")
args <- parser$parse_args()

if (is.null(args[["input"]]) {
    args[["input"]] <- "outputs/fix_nu"
}
source("src/analysis/functions.R")

theme_set(theme_bw())

fit_fix <- readRDS(file.path(args[["input"]], "fix.rds"))
fit_var <- readRDS(file.path(args[["input"]], "var.rds"))

summary_var <- Summary(fit_var)
summary_fix <- Summary(fit_fix)

ord <- order(summary_fix@parameters$mu[, "median"])


g1 <- plot_hpds(summary_var, summary_fix, "mu", ord)
g2 <- plot_hpds(summary_var, summary_fix, "delta", ord)
g3 <- plot_hpds(summary_var, summary_fix, "epsilon", ord)

gg <- ggplotGrob(g1)
legend <- gg$grobs[[grep("guide-box", gg$layout$name)]]

t <- theme(legend.position = "none")

top <- plot_grid(g1 + t, g2 + t, g3 + t, labels = "AUTO", nrow = 1)
combined <- plot_grid(top, legend, nrow = 2, rel_heights = c(0.9, 0.1))
ggsave("figs/fixnu-comparison.pdf", width = 9, height = 4)


d1 <- plot_hpd_diff(summary_var, summary_fix, "mu", ord)
d2 <- plot_hpd_diff(summary_var, summary_fix, "delta", ord)
d3 <- plot_hpd_diff(summary_var, summary_fix, "epsilon", ord)

dop <- plot_grid(d1, d2, d3, labels = "AUTO", nrow = 1)
ggsave("figs/fixnu-comparison-diff.pdf", width = 9, height = 3.7)
