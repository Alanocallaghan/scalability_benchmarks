library("here")
source(here("src/analysis/preamble.R"))

## compare correlation and hpd width
calc_norm_factors <- function(chain) {
    params <- chain@parameters
    if (is.null(params$phi)) {
        params$nu
    } else {
        params$nu * params$phi
    }
}

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
    scale_x_discrete(name = "Number of partitions") +
    scale_y_continuous(
        name = "Pearson correlation",
        limits = c(0.8, 1)
    ) +
    theme(panel.grid = element_blank())

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
    scale_x_discrete(name = "Number of partitions") +
    scale_y_continuous(name = "Pearson correlation", limits = c(0, 1)) +
    theme(panel.grid = element_blank())

ggsave(here("figs/norm_plot_hpd.pdf"), width = 5, height = 5)
