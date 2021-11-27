## comparing scran and BASiCS
library("ggplot2")
library("viridis")
library("scran")
library("ggpointdensity")

fit <- readRDS("outputs/divide_and_conquer/data-chen_nsubsets-1_seed-14_by-gene/chains.rds")
sce <- readRDS("rdata/chen.rds")
sf <- calculateSumFactors(sce)
bf <- colMedians(fit@parameters$nu)

norm_comp_df <- data.frame(scran = sf, BASiCS = bf, data = "chen")

sce_both <- readRDS("rdata/ibarra-soria.rds")
fit <- readRDS("outputs/true-positives/data-ibarra-soria_nsubsets-1_seed-14.rds")
sf <- calculateSumFactors(sce_both[, sce_both$Cell_type == "SM"])
bf <- colMedians(fit$mcmc$sm@parameters$nu)

norm_comp_df <- rbind(norm_comp_df, data.frame(scran = sf, BASiCS = bf, data = "ibarra-som"))

sf <- calculateSumFactors(sce_both[, sce_both$Cell_type == "PSM"])
bf <- colMedians(fit$mcmc$psm@parameters$nu)
norm_comp_df <- rbind(norm_comp_df, data.frame(scran = sf, BASiCS = bf, data = "ibarra-presom"))

r <- range(c(norm_comp_df$scran, norm_comp_df$BASiCS))
# r <- range(c(sf, bf))
# chen_sf_scran <- data.frame(scran = sf, BASiCS = bf)

g <- ggplot(norm_comp_df) +
    aes(scran, BASiCS) +
    geom_pointdensity() +
    facet_wrap(~data) +
    scale_x_log10(limits = r) +
    scale_y_log10(limits = r) +
    scale_colour_viridis(name = "Number of neighbours") +
    labs(x = "scran size factors", y = "BASiCS size factors") +
    geom_text(
        x = 0.3, y = 10,
        hjust = 0,
        mapping = aes(
            label = sprintf(
                "Peason's rho: %0.03f", cor(scran, BASiCS)
            )
        )
    ) +
    theme_bw() +
    theme(
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    theme(legend.position = "bottom")

ggsave("figs/scran_basics.pdf", width = 8, height = 4)


