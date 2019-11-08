get_ess <- function(chain, param) {
  effectiveSize(as.mcmc(chain@parameters[[param]]))
}


plot_all_ess <- function(df, param) {
  ess_all_list <- mclapply(seq_len(nrow(df)), function(i) {
    chain <- readRDS(df[[i, "file"]])
    if (length(chain) > 1) {
      suppressMessages(
        chain <- Scalability:::combine_subposteriors(
          chain,
          subset_by = "gene",
          mc.cores = 1
        )
      )
    }
    data.frame(
      data = df[[i, "data"]],
      chains = df[[i, "chains"]],
      seed = df[[i, "seed"]],
      by = df[[i, "by"]],
      feature = colnames(chain@parameters[[param]]),
      ESS = get_ess(chain, param)
    )
  }, mc.cores = 2)  

  ess_all <- bind_rows(ess_all_list)
  ess_all[which(ess_all[["chains"]] == 1), "by"] <- "Reference"
  ess_all[ess_all[["by"]] == "advi", "by"] <- "ADVI"
  ess_all[ess_all[["by"]] == "gene", "by"] <- "Divide and conquer"
  ess_all[["chains"]] <- factor(ess_all[["chains"]], levels = c(1, 2, 4, 8, 16, 32))

  ess_all[["data"]] <- sub(
    "([[:alpha:]])", "\\U\\1",
    ess_all[["data"]],
    perl = TRUE
  )
  ess_all[["data"]] <- sub(
    "Pbmc",
    "10x PBMC",
    ess_all[["data"]]
  )

  g <- ggplot(ess_all, aes(y = ESS, x = chains, color = by, fill = by)) + 
    # geom_jitter(height = 0, width = 0.2) +
    geom_violin(alpha = 0.2) +
    geom_boxplot(alpha = 0.2, width = 0.1, outlier.colour = NA) +
    facet_wrap(~data, scales = "free_x") +
    scale_fill_brewer(
      name = "Inference method",
      palette = "Dark2",
      aesthetics = c("fill", "color")
    ) +
    labs(
      x = "Partitions",
      y = "Effective sample size"
    )

  ggsave(g, file = paste0("figs/ess_", param, ".pdf"), width = 7, height = 6)
  invisible(g)
}

plot_all_ess(df, "mu")
plot_all_ess(df, "epsilon")
