
stop("Finish the script dickhead!")



calc_norm_factors <- function(chain) {
    params <- chain@parameters
    if (is.null(params$phi)) {
        params$nu
    } else {
        params$nu * params$phi
    }
}

wassersteins <- lapply(seq_len(nrow(df)),
    function(i) {
        cat(i, "/", nrow(df), "\n")
        file <- df[[i, "file"]]
        chain <- readRDS(file)
        match <- which(references$data == df[[i, "data"]])
        ref_norm <- calc_norm_factors(references[[match[[1]], "chain"]])
        my_norm <- calc_norm_factors(chain)
        mean(
            sapply(1:ncol(ref_norm),
            function(j) {
                    transport::wasserstein1d(
                        ref_norm[, j],
                        my_norm[, j]
                    )
                }
            )
        )
    }
)

df$wasserstein <- wassersteins
# sdf <- df %>% filter(!is.na(chains))
# advi_sdf <- sdf %>% filter(is.na(chains)) 


g <- ggplot(df,
  aes(
    x = factor(chains),
    y = wasserstein
  )
) +
  geom_quasirandom(
    groupOnX = TRUE,
    dodge.width = 0.5,
    # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
    size = 0.5
  ) +
  # geom_point(
  #   position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
  #   size = 0.5
  # ) +
  # geom_violin() +
  facet_wrap(~data, nrow = 2, ncol = 2) +
  scale_x_discrete(name = "Partitions") +
  scale_y_continuous(name = "Portion of genes perturbed", labels = scales::percent) +
  theme(text = element_text(size = 18))

ggsave(here("figs/diffexp_plot.pdf"), width = 12, height = 8)
