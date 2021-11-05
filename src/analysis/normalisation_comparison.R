
## compare correlation and hpd width
calc_norm_factors <- function(chain) {
    params <- chain@parameters
    if (is.null(params$phi)) {
        params$nu
    } else {
        params$nu * params$phi
    }
}

# dr <- df[df$data == "chen" & df$chains == 1, ]
# dd <- df[df$data == "chen" & df$chains == 4, ]
# x1 <- readRDS(dr$file[[1]])
# x2 <- readRDS(dd$file[[1]])
# cm1 <- colMedians(calc_norm_factors(x1))
# cm2 <- colMedians(calc_norm_factors(x2))
# # plot(cm1, cm2)

cors <- lapply(seq_len(nrow(df)),
    function(i) {
        cat(i, "/", nrow(df), "\n")
        file <- df[[i, "file"]]
        chain <- readRDS(file)
        match <- which(references$data == df[[i, "data"]])
        ref_norm <- calc_norm_factors(references[[match[[1]], "chain"]])
        my_norm <- calc_norm_factors(chain)
        med_r <- colMedians(ref_norm)
        my_r <- colMedians(my_norm)
        ref_hpd <- coda::HPDinterval(as.mcmc(ref_norm))
        my_hpd <- coda::HPDinterval(as.mcmc(my_norm))
        data.frame(
            cor_mean = cor(med_r, my_r),
            cor_lower = cor(ref_hpd[, "lower"], my_hpd[, "lower"]),
            cor_upper = cor(ref_hpd[, "upper"], my_hpd[, "upper"])
        )
    }
)
df_norm <- df
cor_df <- do.call(rbind, cors)
df_norm <- cbind(df_norm, cor_df)
df_norm <- df_norm[df_norm$chains != 1, ]

# sdf <- df %>% filter(!is.na(chains))
# advi_sdf <- sdf %>% filter(is.na(chains)) 

g <- ggplot(
        df_norm,
        aes(
            x = factor(chains),
            # colour = by,
            y = cor_mean
        )
    ) +
    geom_quasirandom(
        groupOnX = TRUE,
        dodge.width = 0.5,
        # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
        size = 0.5
    ) +
    # geom_point(
    #     position = position_jitterdodge(
    #         jitter.width = 0.1,
    #         jitter.height = 0
    #     ),
    #     size = 0.5
    # ) +
    # geom_violin() +
    facet_wrap(~data, nrow = 2, ncol = 2) +
    scale_x_discrete(name = "Partitions") +
    scale_y_continuous(
    name = "Pearson correlation"
    ) +
    theme(text = element_text(size = 18))

ggsave(here("figs/norm_plot.pdf"), width = 5, height = 5)


df_norm
g <- ggplot(
        df_norm,
        aes(
            x = factor(chains),
        )
    ) +
    geom_quasirandom(
        mapping = aes(colour = "95% HPD interval (lower)", y = cor_lower),
        groupOnX = TRUE,
        dodge.width = 0.5,
        # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
        size = 0.5
    ) +
    geom_quasirandom(
        mapping = aes(colour = "95% HPD interval (upper)", y = cor_upper),
        groupOnX = TRUE,
        dodge.width = 0.5,
        # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
        size = 0.5
    ) +
    # geom_point(
    #     position = position_jitterdodge(
    #         jitter.width = 0.1,
    #         jitter.height = 0
    #     ),
    #     size = 0.5
    # ) +
    # geom_violin() +
    facet_wrap(~data, nrow = 2, ncol = 2) +
    scale_x_discrete(name = "Partitions") +
    scale_y_continuous(name = "Pearson correlation") +
    theme(text = element_text(size = 18))

ggsave(here("figs/norm_plot_hpd.pdf"), width = 5, height = 5)



## comparing scran and BASiCS
library("ggplot2")
library("viridis")
library("ggpointdensity")
fit <- readRDS("outputs/divide_and_conquer/data-chen_nsubsets-1_seed-14_by-gene/chains.rds")
sce <- readRDS("rdata/chen.rds")
sf <- calculateSumFactors(sce)
bf <- colMedians(fit@parameters$nu)
r <- range(c(sf, bf))
chen_sf_scran <- data.frame(scran = sf, BASiCS = bf)
g <- ggplot(chen_sf_scran) +
    aes(scran, BASiCS) +
    # geom_point() +
    geom_pointdensity() +
    scale_x_log10(limits = r) +
    scale_y_log10(limits = r) +
    scale_colour_viridis(name = "Number of neighbours") +
    labs(x = "scran size factors", y = "BASiCS size factors") +
    annotate("text", x = 3.1, y = 1.2,
        hjust = 0,
        label = sprintf(
            "Peason's rho: %0.03f", cor(chen_sf_scran[[1]], chen_sf_scran[[2]])
        )
    ) +
    theme_bw() +
    theme(
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    theme(legend.position = "bottom")

ggsave("figs/chen_scran_basics.pdf", width = 4.5, height = 4)



## wasserstein idea
# wassersteins <- lapply(seq_len(nrow(df)),
#     function(i) {
#         cat(i, "/", nrow(df), "\n")
#         file <- df[[i, "file"]]
#         chain <- readRDS(file)
#         match <- which(references$data == df[[i, "data"]])
#         ref_norm <- calc_norm_factors(references[[match[[1]], "chain"]])
#         my_norm <- calc_norm_factors(chain)
#         mean(
#             sapply(1:ncol(ref_norm),
#             function(j) {
#                     transport::wasserstein1d(
#                         ref_norm[, j],
#                         my_norm[, j]
#                     )
#                 }
#             )
#         )
#     }
# )

# df$wasserstein <- wassersteins
# # sdf <- df %>% filter(!is.na(chains))
# # advi_sdf <- sdf %>% filter(is.na(chains)) 


# g <- ggplot(df,
#   aes(
#     x = factor(chains),
#     y = wasserstein
#   )
# ) +
#   geom_quasirandom(
#     groupOnX = TRUE,
#     dodge.width = 0.5,
#     # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
#     size = 0.5
#   ) +
#   # geom_point(
#   #   position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
#   #   size = 0.5
#   # ) +
#   # geom_violin() +
#   facet_wrap(~data, nrow = 2, ncol = 2) +
#   scale_x_discrete(name = "Partitions") +
#   scale_y_continuous(name = "Portion of genes perturbed", labels = scales::percent) +
#   theme(text = element_text(size = 18))

# ggsave(here("figs/norm_plot.pdf"), width = 12, height = 8)

