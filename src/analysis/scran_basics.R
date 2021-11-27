## comparing scran and BASiCS
library("ggplot2")
library("cowplot")
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


gs <- lapply(
    unique(norm_comp_df$data),
    function(d) {
        norm_comp_sub <- norm_comp_df[norm_comp_df$data == d, ]
        ggplot(norm_comp_sub) +
            aes(scran, BASiCS) +
            geom_pointdensity() +
            scale_x_log10(limits = r) +
            scale_y_log10(limits = r) +
            scale_colour_viridis(name = "Number of neighbours") +
            labs(x = "scran size factors", y = "BASiCS size factors") +
            annotate(
                geom = "text",
                x = 0.3, y = 10.1,
                hjust = 0,
                vjust = 0,
                label = sprintf(
                    "Peason's rho: %0.03f",
                    cor(norm_comp_sub$scran, norm_comp_sub$BASiCS)
                )
            ) +
            theme_bw() +
            theme(
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            theme(legend.position = "none")
    }
)
plot_grid(plotlist = gs, nrow = 1, labels = "AUTO")
ggsave("figs/scran_basics.pdf", width = 8, height = 3.5)
