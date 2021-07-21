get_normalised_hpd_width <- function(chain, param) {
  hpd <- HPDinterval(as.mcmc(chain@parameters[[param]]))
  (hpd[, "upper"] - hpd[, "lower"]) / rowMeans(hpd)
}

plot_all_hpds <- function(df, param) {
  hpds_all <- mclapply(
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
        hpd = get_normalised_hpd_width(chain, param)
      )
    }, mc.cores = 2
  )  

  hpdf_all <- bind_rows(hpds_all)
  hpdf_all[which(hpdf_all[["chains"]] == 1), "by"] <- "Reference"
  hpdf_all[hpdf_all[["by"]] == "advi", "by"] <- "ADVI"
  hpdf_all[hpdf_all[["by"]] == "gene", "by"] <- "Divide and conquer"
  hpdf_all[["chains"]] <- factor(hpdf_all[["chains"]], levels = c(1, 2, 4, 8, 16, 32))


  hpdf_ordered <- hpdf_all %>% 
    group_by(data, chains, seed, by) %>% 
    arrange(feature, .by_group = TRUE)

  hpdf_ref <- hpdf_ordered[which(hpdf_ordered$chains == 1 & !is.na(hpdf_ordered$chains)), ]
  hpdf_nr <- hpdf_ordered[which(hpdf_ordered$chains != 1 | is.na(hpdf_ordered$chains)), ]

  hpdf_nr <- as.data.frame(hpdf_nr)
  hpdf_ref <- as.data.frame(hpdf_ref)

  for (dataset in unique(hpdf_ref[["data"]])) {
    ind_nr <- hpdf_nr[["data"]] == dataset
    ind_r <- hpdf_ref[["data"]] == dataset
    d1 <- hpdf_nr[ind_nr, ]
    d2 <- hpdf_ref[ind_r, ]
    stopifnot(all(d1$data == d2$data))
    stopifnot(all(d1$feature == d2$feature))
    # stopifnot(all(hpdf_nr[ind_nr, "data"] == hpdf_ref[ind_r, "data"]))
    # stopifnot(all(hpdf_nr[ind_nr, "feature"] == hpdf_ref[ind_r, "feature"]))
    hpdf_nr[ind_nr, "hpd"] <- hpdf_nr[ind_nr, "hpd"] - hpdf_ref[ind_r, "hpd"]
  }


  hpdf_nr[["data"]] <- sub(
    "([[:alpha:]])", "\\U\\1",
    hpdf_nr[["data"]],
    perl = TRUE
  )
  hpdf_nr[["data"]] <- sub(
    "Pbmc",
    "10x PBMC",
    hpdf_nr[["data"]]
  )


  g <- ggplot(hpdf_nr, aes(y = hpd, x = chains, color = by, fill = by)) + 
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
      y = "Normalised HPD width relative to AMWG normalised HPD width"
    )

  ggsave(g, file = paste0("figs/hpd_width_", param, ".pdf"), width = 7, height = 6)
  invisible(g)
}
source(here("src/analysis/preamble.R"))

plot_all_hpds(df, "mu")
plot_all_hpds(df, "epsilon")
