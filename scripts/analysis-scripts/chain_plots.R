do_fit_plot <- function(i, maxdf, references) {
  c <- readRDS(maxdf[[i, "file"]])
  if (is.list(c)) {  
    suppressMessages(
      c <- combine_subposteriors(
        c,
        subset_by = "gene"
      )
    )
  }
  d <- maxdf[[i, "data"]]
  rc <- references[[which(references$data == d), "chain"]]
  c <- offset_correct(rc, c)
  g1 <- BASiCS_ShowFit(rc) + ggtitle("Reference")
  g2 <- BASiCS_ShowFit(c) + ggtitle("Divide and conquer")
  cowplot::plot_grid(g1, g2)
}

do_de_plot <- function(i, maxdf, references) {
  d <- maxdf[[i, "data"]]
  rc <- references[[which(references$data == d), "chain"]]

  c <- readRDS(maxdf[[i, "file"]])
  if (is.list(c)) {  
    suppressMessages(
      c <- combine_subposteriors(
        c,
        gene_order = rownames(rc),
        cell_order = colnames(rc),
        subset_by = "gene"
      )
    )
  }
  de <- BASiCS_TestDE(rc, c,
    GroupLabel1 = "Ref",
    GroupLabel2 = "D & C",
    Plot = FALSE,
    PlotOffset = FALSE,
    EFDR_M = NULL,
    EFDR_D = NULL,
    EFDR_R = NULL
  )
  BASiCS_PlotDE(de, Which = c("MAPlot", "VolcanoPlot"))
}





## Mean-variance curves and DE plots for worst for each data
maxdfm <- df %>% filter(!is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)


maxdfe <- df %>% filter(!is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)


fit_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_fit_plot,
  maxdf = maxdfe,
  references = references
)

## MA plot of worst case for worst for each data
de_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_de_plot,
  maxdf = maxdfe,
  references = references
)


## Posterior variance plots
maxdfadvim <- df %>% filter(is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)

maxdfadvie <- df %>% filter(is.na(chains)) %>% 
  group_by(data) %>% 
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)


fit_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_fit_plot,
  maxdf = maxdfadvie, 
  references = references
)

## MA plot of worst case for worst for each data
de_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_de_plot,
  maxdf = maxdfadvie,
  references = references
)




hpd_plots <- lapply(seq_len(nrow(maxdfadvie)),
  function(j) {
    c <- readRDS(maxdfadvie[[j, "file"]])
    d <- maxdf[[j, "data"]]
    rc <- references[[which(references$data == d), "chain"]]
    c <- offset_correct(rc, c)
    list(
      plot_hpd_interval(
        rc,
        "Reference",
        c,
        "ADVI",
        "mu",
        type = "ma",
        log = FALSE,
        bins = 50
      ) + 
        labs(
          title = NULL,
          x = bquote(mu[i]),
          y = "log2(ADVI HPD interval width / Reference HPD interval width)"
        ) +
        geom_hline(aes(yintercept = 0)),
      plot_hpd_interval(
        rc,
        "Reference",
        c,
        "ADVI",
        "epsilon",
        type = "ma",
        log = FALSE,
        bins = 50
      ) +
        labs(
          title = NULL,
          x = bquote(mu[i]),
          y = "log2(ADVI HPD interval width / Reference HPD interval width)"
        ) +
        geom_hline(aes(yintercept = 0))
    )
  }
)
