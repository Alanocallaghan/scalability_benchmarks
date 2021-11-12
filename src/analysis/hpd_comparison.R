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
parser$add_argument("-d", "--data")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

data <- args[["dataset"]]
if (data == "ibarra-soria") {
    args <- list(
        reference = "outputs/true-positives/data-ibarra-soria_nsubsets-1_seed-14.rds",
        divide = "outputs/true-positives/data-ibarra-soria_nsubsets-16_seed-14.rds",
        advi = "outputs/true-positives/data-ibarra-soria_advi-14.rds",
        output =  "figs/hpd/ibarra-soria"
    )
} else {
    args <- list(
        reference = sprintf("outputs/divide_and_conquer/data-%s_nsubsets-16_seed-14_by-gene/chains.rds", data),
        divide = sprintf("outputs/divide_and_conquer/data-%s_nsubsets-1_seed-14_by-gene/chains.rds", data),
        advi = sprintf("outputs/advi/data-%s_seed-14/chain.rds", data),
        output = sprintf("figs/hpd/%s", data)
    )
}

w <- 8
h <- 3

output <- args[["output"]]
dir.create(output, showWarnings = FALSE)

source("src/analysis/functions.R")


theme_set(theme_bw())

fit_ref <- readRDS(args[["reference"]])
fit_dc <- readRDS(args[["divide"]])
fit_advi <- readRDS(args[["advi"]])

if (data == "ibarra-soria") {
    fit_ref <- fit_ref$mcmc[[1]]
    fit_dc <- fit_dc$mcmc[[1]]
    fit_advi <- fit_advi$mcmc[[1]]
}

summary_ref <- Summary(fit_ref)
summary_dc <- Summary(fit_dc)
summary_advi <- Summary(fit_advi)

ord <- order(summary_ref@parameters$mu[, "median"])

################################################################################
## HPD and point estimates across inference methods
################################################################################

ord <- order(summary_ref@parameters$mu[, "median"])

g1 <- plot_hpds(summary_ref, summary_dc, "mu", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "Divide and conquer")
)
g2 <- plot_hpds(summary_ref, summary_dc, "delta", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "Divide and conquer")
)
g3 <- plot_hpds(summary_ref, summary_dc, "epsilon", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "Divide and conquer")
)

gg <- ggplotGrob(g1)
legend <- gg$grobs[[grep("guide-box", gg$layout$name)]]

t <- theme(legend.position = "none")

top <- plot_grid(g1 + t, g2 + t, g3 + t, labels = "AUTO", nrow = 1)
combined <- plot_grid(top, legend, nrow = 2, rel_heights = c(0.9, 0.1))
ggsave(file.path(output, "divide-and-conquer-hpd.pdf"), width = w, height = h)



g1 <- plot_hpds(summary_ref, summary_advi, "mu", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "ADVI")
)
g2 <- plot_hpds(summary_ref, summary_advi, "delta", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "ADVI")
)
g3 <- plot_hpds(summary_ref, summary_advi, "epsilon", ord,
    scalename = "Inference method     ",
    labels = c("Reference", "ADVI")
)

gg <- ggplotGrob(g1)
legend <- gg$grobs[[grep("guide-box", gg$layout$name)]]

t <- theme(legend.position = "none")

top <- plot_grid(g1 + t, g2 + t, g3 + t, labels = "AUTO", nrow = 1)
combined <- plot_grid(top, legend, nrow = 2, rel_heights = c(0.9, 0.1))

ggsave(file.path(output, "advi-hpd.pdf"), width = w, height = h)




################################################################################
## Difference in HPDs between inference methods
################################################################################



d1 <- plot_hpd_diff(summary_ref, summary_dc, "mu", ord)
d2 <- plot_hpd_diff(summary_ref, summary_dc, "delta", ord)
d3 <- plot_hpd_diff(summary_ref, summary_dc, "epsilon", ord)

dop <- plot_grid(d1, d2, d3, labels = "AUTO", nrow = 1)

ggsave(file.path(output, "divide-and-conquer-hpd-diff.pdf"), width = w, height = h)


d1 <- plot_hpd_diff(summary_ref, summary_advi, "mu", ord)
d2 <- plot_hpd_diff(summary_ref, summary_advi, "delta", ord)
d3 <- plot_hpd_diff(summary_ref, summary_advi, "epsilon", ord)

dop <- plot_grid(d1, d2, d3, labels = "AUTO", nrow = 1)


ggsave(file.path(output, "advi-hpd-diff.pdf"), width = w, height = h)
