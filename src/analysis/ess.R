library("here")

get_ess <- function(chain, param) {
  effectiveSize(as.mcmc(chain@parameters[[param]]))
}

plot_all_ess <- function(df, param) {
  ess_all_list <- mclapply(
    seq_len(nrow(df)),
    function(i) {
      cat(i, "/", nrow(df), "\n")
      chain <- readRDS(df[[i, "file"]])
      if (length(chain) > 1) {
        suppressMessages(
          chain <- BASiCS:::.combine_subposteriors(
            chain,
            SubsetBy = "gene"
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
    },
    mc.cores = 5
  )

  ess_all <- bind_rows(ess_all_list)
  ess_all[which(ess_all[["chains"]] == 1), "by"] <- "Reference"
  ess_all <- ess_all[ess_all[["by"]] != "advi", ]
  ess_all[ess_all[["by"]] == "gene", "by"] <- "Divide and conquer"
  ess_all[["chains"]] <- factor(
    ess_all[["chains"]],
    levels = sort(unique(ess_all[["chains"]]))
  )

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

  g <- ggplot(ess_all, aes(y = ESS, x = chains
      # , color = by, fill = by
      )
    ) +
    # geom_jitter(height = 0, width = 0.2) +
    geom_violin(alpha = 0.2) +
    geom_boxplot(alpha = 0.2, width = 0.1, outlier.colour = NA) +
    facet_wrap(~data, scales = "free_x") +
    scale_y_log10() +
    # scale_fill_brewer(
    #   name = "Inference method",
    #   palette = "Dark2",
    #   aesthetics = c("fill", "color")
    # ) +
    labs(
      x = "Partitions",
      y = "Effective sample size"
    )

  ggsave(g, file = sprintf("figs/ess_%s.pdf", param), width = 7, height = 6)
  invisible(g)
}
source(here("src/analysis/preamble.R"))

plot_all_ess(df, "mu")
plot_all_ess(df, "epsilon")
