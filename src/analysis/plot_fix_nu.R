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

args[["input"]] <- "outputs/fix_nu"

theme_set(theme_bw())

fit_fix <- readRDS(file.path(args[["input"]], "fix.rds"))
fit_var <- readRDS(file.path(args[["input"]], "var.rds"))

summary_var <- Summary(fit_var)
summary_fix <- Summary(fit_fix)

ord <- order(summary_fix@parameters$mu[, "median"])

plot_hpds <- function(summary_var, summary_fix, param = "mu", ord) {
    df_var <- as.data.frame(summary_var@parameters[[param]])
    df_fix <- as.data.frame(summary_fix@parameters[[param]])
    # ord <- order(df_var$median)
    df_var <- df_var[ord, ]
    df_fix <- df_fix[ord, ]
    df_var$index <- 1:nrow(df_var)
    df_fix$index <- 1:nrow(df_fix)

    if (param %in% c("mu", "delta")) {
        scale <- scale_y_log10(name = param)
    } else {
        scale <- scale_y_continuous(name = param)
    }

    ggplot() +
        geom_line(
            data = df_var,
            alpha = 0.6,
            aes(x = index, y = median, colour = "Inferred")
        ) +
        geom_line(
            data = df_fix,
            alpha = 0.6,
            aes(x = index, y = median, colour = "Fixed")
        ) +
        geom_ribbon(
            data = df_var,
            alpha = 0.4,
            aes(x = index, ymin = lower, ymax = upper, fill = "Inferred")
        ) +
        geom_ribbon(
            data = df_fix,
            alpha = 0.4,
            aes(x = index, ymin = lower, ymax = upper, fill = "Fixed")
        ) +
        scale +
        xlab("Gene") +
        guides(colour = "none") +
        theme(legend.position = "bottom") +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        ) +
        scale_colour_manual(
            name = "Normalisation    ",
            values = c("firebrick", "dodgerblue"),
            aesthetics = c("fill", "colour")
        )
}

g1 <- plot_hpds(summary_var, summary_fix, "mu", ord)
g2 <- plot_hpds(summary_var, summary_fix, "delta", ord)
g3 <- plot_hpds(summary_var, summary_fix, "epsilon", ord)

gg <- ggplotGrob(g1)
legend <- gg$grobs[[grep("guide-box", gg$layout$name)]]

t <- theme(legend.position = "none")

top <- plot_grid(g1 + t, g2 + t, g3 + t, labels = "AUTO", nrow = 1)
combined <- plot_grid(top, legend, nrow = 2, rel_heights = c(0.9, 0.1))
ggsave("figs/fixnu-comparison.pdf", width = 9, height = 4)


plot_hpd_diff <- function(summary_var, summary_fix, param = "mu", ord) {
    df_var <- as.data.frame(summary_var@parameters[[param]])
    df_fix <- as.data.frame(summary_fix@parameters[[param]])
    # ord <- order(df_var$median)
    df_var <- df_var[ord, ]
    df_fix <- df_fix[ord, ]
    df_diff <- df_var - df_fix
    df_diff$index <- 1:nrow(df_var)

    ggplot() +
        geom_ribbon(
            data = df_diff,
            colour = "darkgoldenrod",
            alpha = 0.5,
            aes(x = index, ymin = lower, ymax = upper)
        ) +
        geom_line(
            data = df_diff,
            alpha = 0.8,
            aes(x = index, y = median)
        ) +
        xlab("Gene") +
        ylab(sprintf("Difference in posterior quantities (%s)", param)) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

        # +
        # guides(colour = "none") +
        # theme(legend.position = "bottom") +
        # scale_colour_manual(
        #     name = "Normalisation    ",
        #     values = c("firebrick", "dodgerblue"),
        #     aesthetics = c("fill", "colour")
        # )
}

d1 <- plot_hpd_diff(summary_var, summary_fix, "mu", ord)
d2 <- plot_hpd_diff(summary_var, summary_fix, "delta", ord)
d3 <- plot_hpd_diff(summary_var, summary_fix, "epsilon", ord)

dop <- plot_grid(d1, d2, d3, labels = "AUTO", nrow = 1)
ggsave("figs/fixnu-comparison-diff.pdf", width = 9, height = 3.5)
