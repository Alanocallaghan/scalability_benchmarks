
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
        var_r <- colVars(ref_norm)
        var_my <- colVars(my_norm)
        ref_hpd <- coda::HPDinterval(as.mcmc(ref_norm))
        my_hpd <- coda::HPDinterval(as.mcmc(my_norm))
        data.frame(
            cor_mean = cor(med_r, my_r),
            cor_var = cor(var_my, var_r),
            cor_lower = cor(ref_hpd[, "lower"], my_hpd[, "lower"]),
            cor_upper = cor(ref_hpd[, "upper"], my_hpd[, "upper"])
        )
    }
)
df_norm <- df
cor_df <- do.call(rbind, cors)
df_norm <- cbind(df_norm, cor_df)
df_norm$data <- gsub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    df_norm$data,
    perl = TRUE
)
df_norm_advi <- df_norm[df_norm$by == "advi", ]
df_norm <- df_norm[which(df_norm$chains != 1), ]
df_norm_advi_sum <- df_norm_advi %>%
    group_by(data) %>%
    summarise(
        mean_cor_mean = mean(cor_mean),
        mean_cor_var = mean(cor_var)
    )
df_norm_advi <- df_norm_advi %>%
    group_by(data) %>%
    filter(row_number() == 1) %>%
    cbind(df_norm_advi_sum[-1])
df_norm_advi <- df_norm_advi[df_norm_advi$data != "Chen", ]
# sdf <- df %>% filter(!is.na(chains))
# advi_sdf <- sdf %>% filter(is.na(chains)) 

g <- ggplot() +
    geom_quasirandom(
        data = df_norm,
        mapping = aes(
            x = factor(chains),
            # colour = by,
            y = cor_mean
        ),
        groupOnX = TRUE,
        dodge.width = 0.5,
        size = 0.5
    ) +
    geom_hline(
        data = df_norm_advi,
        aes(yintercept = mean_cor_mean),
        linetype = "dashed"
    ) +
    facet_wrap(~data, nrow = 2, ncol = 2) +
    scale_x_discrete(name = "Partitions") +
    scale_y_continuous(
        name = "Pearson correlation",
        limits = c(0.8, 1)
    ) +
    theme(text = element_text(size = 18), panel.grid = element_blank())

ggsave(here("figs/norm_plot.pdf"), width = 5, height = 5)

g <- ggplot(
        df_norm,
        aes(
            x = factor(chains),
        )
    ) +
    geom_quasirandom(
        mapping = aes(y = cor_var),
        groupOnX = TRUE,
        dodge.width = 0.5,
        size = 0.5
    ) +
    geom_hline(
        data = df_norm_advi,
        aes(yintercept = mean_cor_var),
        linetype = "dashed"
    ) +
    facet_wrap(~data, nrow = 2, ncol = 2) +
    scale_x_discrete(name = "Partitions") +
    scale_y_continuous(name = "Pearson correlation", limits = c(0, 1)) +
    theme(text = element_text(size = 18), panel.grid = element_blank())

ggsave(here("figs/norm_plot_hpd.pdf"), width = 5, height = 5)


## comparing scran and BASiCS
library("ggplot2")
library("viridis")
library("scran")
library("ggpointdensity")

fit <- readRDS("outputs/divide_and_conquer/data-chen_nsubsets-1_seed-14_by-gene/chains.rds")
sce <- readRDS("rdata/chen.rds")
sf <- calculateSumFactors(sce)
bf <- colMedians(fit@parameters$nu)

df <- data.frame(scran = sf, basics = bf, data = "chen")

sce_both <- readRDS("rdata/ibarra-soria.rds")
fit <- readRDS("outputs/true-positives/data-ibarra-soria_nsubsets-1_seed-14.rds")
sf <- calculateSumFactors(sce_both[, sce_both$Cell_type == "SM"])
bf <- colMedians(fit$mcmc$sm@parameters$nu)

df <- rbind(df, data.frame(scran = sf, basics = bf, data = "ibarra-som"))

sf <- calculateSumFactors(sce_both[, sce_both$Cell_type == "PSM"])
bf <- colMedians(fit$mcmc$psm@parameters$nu)
df <- rbind(df, data.frame(scran = sf, basics = bf, data = "ibarra-presom"))

r <- range(c(norm_df$scran, norm_df$basics))
# r <- range(c(sf, bf))
# chen_sf_scran <- data.frame(scran = sf, BASiCS = bf)
g <- ggplot(chen_sf_scran) +
    aes(scran, BASiCS) +
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
#   theme(text = element_text(size = 18))

# ggsave(here("figs/norm_plot.pdf"), width = 12, height = 8)

