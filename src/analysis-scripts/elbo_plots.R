dir.create("figs/elbo", recursive = TRUE, showWarnings = FALSE)

parsed_elbos <- lapply(advi_elbo, parse_elbo)


elbo_df <- as.data.frame(advi_df)
elbo_df$elbos <- parsed_elbos
elbo_df$elbos <- lapply(seq_len(nrow(elbo_df)), 
  function(i) {
    d <- elbo_df$elbos[[i]]
    d$data <- elbo_df[[i, "data"]]
    d$seed <- elbo_df[[i, "seed"]]
    d
  }
)

elbo_df <- do.call(rbind, elbo_df$elbos)
elbo_df$Data <- sub(
  "([[:alpha:]])", "\\U\\1",
  elbo_df$data,
  perl = TRUE
)
elbo_df$Data <- sub(
  "Pbmc",
  "10x PBMC",
  elbo_df$Data
)
plots <- lapply(unique(elbo_df$data), 
  function(d) {
    df <- elbo_df[elbo_df$data == d, ]
    D <- df$Data[[1]]
    g <- ggplot(df, aes(x = iter, y = abs(ELBO), color = factor(seed))) +
      geom_line() +
      theme(legend.position = "none") +
      scale_y_log10() +
      # ggtitle(D) +
      labs(x = "Iteration")
    ggsave(g, file = here(paste0("figs/elbo/", d, ".pdf")), width = 7, height = 7)
    g
  }
)
