read_triplets <- function(triplets, combine = FALSE) {
  rows <- lapply(
    triplets,
    function(x) {
      row <- readRDS(x[[2]])
      row <- as.data.frame(row)
      row[["time"]] <- readRDS(x[[3]])[["elapsed"]]
      row
    }
  )
  df <- do.call(rbind, rows)
  df[["file"]] <- lapply(triplets, function(x) x[[1]])
  df
}

file2triplets <- function(files) {
  lapply(files, list.files, full.names = TRUE)
}

## Convergence first
parse_elbo <- function(c) {
  c <- gsub("Chain 1:\\s+", "", c)
  normal <- grep("Drawing", c)
  abnormal <- grep("Informational", c)
  if (!length(normal)) {
    stop("Something is deeply wrong here")
  }
  ind_end <- normal - 2
  if (length(abnormal)) {
    ind_end <- abnormal - 1
  }
  elbo <- c[(grep("Begin stochastic", c) + 1):ind_end]
  ## This ain't quite right
  elbo <- gsub("MAY BE DIVERGING... INSPECT ELBO", "", elbo, fixed = TRUE)
  elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
    elbo[-c(1, grep("CONVERGED", elbo))],
    "NOTCONVERGED"
  )
  ## If both mean and median elbo converge this ends up with CONVERGED CONVERGED
  ## and therefore another column
  elbo <- gsub("(MEDIAN |MEAN )?ELBO CONVERGED", "CONVERGED", elbo)
  elbo <- gsub("(\\s+CONVERGED){2}", " CONVERGED", elbo)
  elbo <- strsplit(elbo, "\\s+")
  elbo <- do.call(rbind, elbo)
  colnames(elbo) <- elbo[1, ]
  elbo <- elbo[-1, ]
  elbo <- as.data.frame(elbo, stringsAsFactors = FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}

do_de <- function(
    df,
    ref_df,
    match_column,
    data_dims,
    mc.cores = options("mc.cores")
  ) {

  edr <- lapply(
    seq_len(nrow(df)),
    function(i) {
      cat(i, "/", nrow(df), "\n")
      if (isTRUE(df[[i, "chains"]] == 1)) {
        return(rep(list(NULL), 3))
      }
      ind <- ref_df[[match_column]] == df[[i, match_column]] &
        ref_df[["data"]] == df[i, "data"]
      chain <- readRDS(df[[i, "file"]])
      if (length(chain) > 1) {
        suppressMessages(
          chain <- BASiCS:::.combine_subposteriors(
            chain,
            SubsetBy = "gene"
          )
        )
      }
      cp <- intersect(c("nu", "s", "phi"), names(chain@parameters))
      chain@parameters[cp] <- lapply(
        chain@parameters[cp],
        function(x) {
          colnames(x) <- gsub("_Batch.*", "", colnames(x))
          x
        }
      )
      chain@parameters <- BASiCS:::.reorder_params(
        chain@parameters,
        GeneOrder = rownames(ref_df[[which(ind)[[1]], "chain"]])
      )
      nsamples <- nrow(
        ref_df[[which(ind)[[1]], "chain"]]@parameters[["mu"]]
      )
      params <- setdiff(names(chain@parameters), "RBFLocations")
      chain@parameters[params] <- lapply(
        chain@parameters[params],
        function(x) {
          x[seq_len(nsamples), ]
        }
      )
      suppressMessages(
        de <- BASiCS_TestDE(
          ref_df[[which(ind)[[1]], "chain"]],
          chain,
          Plot = FALSE,
          PlotOffset = FALSE,
          EFDR_M = NULL,
          EFDR_D = NULL,
          EFDR_R = NULL
        )
      )
      lapply(
        de@Results,
        function(x) {
          l <- as.data.frame(x)[["GeneName"]]
          if (!length(l)) NULL else l
        }
      )
    }
    # ,
    # mc.cores = mc.cores
  )

  edr_df <- do.call(rbind, edr)
  colnames(edr_df) <- c("DiffExp", "DiffDisp", "DiffResDisp")
  edr_df <- as.data.frame(edr_df)
  df <- cbind(df, edr_df)

  df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
    df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
    function(x) sapply(x, length)
  )

  df <- merge(df, data_dims)
  df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
    df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / df[["nGenes"]],
    digits = 3
  )
  df
}

jaccard <- function(a, b) {
  a <- na.omit(a)
  b <- na.omit(b)
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  if (union == 0) 0 else intersection / union
}

sourceme <- function(x) {
  cat("Running", x, "\n")
  source(x)
}

do_fit_plot <- function(i, maxdf, references) {
  cat(i, "/", nrow(maxdf), "\n")
  maxdf <- maxdf[i, ]
  c <- readRDS(maxdf[["file"]][[1]])
  b <- maxdf[["by"]]
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        SubsetBy = "gene"
      )
    )
  }
  d <- maxdf[["data"]]
  rc <- references[[which(references$data == d)[[1]], "chain"]]
  l2 <- if (b == "advi") "ADVI" else "Divide and conquer"
  c <- BASiCS:::.offset_correct(rc, c)
  g1 <- BASiCS_ShowFit(rc) 
  g2 <- BASiCS_ShowFit(c) 
  ggsave(g1, 
    file = here(sprintf("figs/fit/%s_ref.pdf", d)),
    width = 4,
    height = 3
  )
  ggsave(g2,
    file = here(sprintf("figs/fit/%s_%s.pdf", d, b)),
    width = 4,
    height = 3
  )
  cowplot::plot_grid(
    g1 + ggtitle("Reference"),
    g2 + ggtitle(l2)
  )
}

do_de_plot <- function(i, maxdf, references) {
  cat(i, "/", nrow(maxdf), "\n")
  maxdf <- maxdf[i, ]
  d <- maxdf[["data"]]
  b <- maxdf[["by"]]
  rc <- references[[which(references$data == d)[[1]], "chain"]]
  cc <- readRDS(maxdf[["file"]][[1]])
  if (is.list(cc)) {
    suppressMessages(
      cc <- BASiCS:::.combine_subposteriors(
        cc,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
    )
  }
  cc@parameters <- BASiCS:::.reorder_params(
    cc@parameters,
    GeneOrder = rownames(rc)
  )

  l2 <- if (b == "advi") "ADVI" else "D & C"
  suppressMessages(
    de <- BASiCS_TestDE(
      rc, cc,
      GroupLabel1 = "Ref",
      GroupLabel2 = l2,
      Plot = FALSE,
      PlotOffset = FALSE,
      EFDR_M = NULL,
      EFDR_D = NULL,
      EFDR_R = NULL
    )
  )
  g1 <- BASiCS_PlotDE(de@Results[[1]],
    Plots = "MA",
    Mu = de@Results$Mean@Table$MeanOverall
  )
  g2 <- BASiCS_PlotDE(de@Results[[1]],
    Plots = "Volcano",
    Mu = de@Results$Mean@Table$MeanOverall
  ) + theme(legend.position = "bottom")
  gg2 <- ggplotGrob(g2)
  legend <- gg$grobs[[grep("guide-box", gg$layout$name)]]
  t <- theme(legend.position = "none")
  together <- cowplot::plot_grid(g1 + t, g2 + t, labels = "AUTO")
  all <- cowplot::plot_grid(together, legend, nrow = 2, rel_height = 0.9, 0.1)
  ggsave(all,
    file = here(sprintf("figs/de/mu_%s_%s.pdf", d, b)),
    width = 6, height = 5
  )
  g <- BASiCS_PlotDE(
    de@Results[[3]],
    Plots = c("MA", "Volcano"),
    Mu = de@Results$Mean@Table$MeanOverall
  )
  ggsave(g,
    file = here(sprintf("figs/de/epsilon_%s_%s.pdf", d, b)),
    width = 6, height = 5
  )
  g
}

plot_hpd_interval <- function(
    chain1,
    xname,
    chain2,
    yname,
    param,
    type = c("ma", "std"),
    log = FALSE,
    ...
  ) {

  x <- extract_hpd_interval(chain1, param)
  y <- extract_hpd_interval(chain2, param)

  df <- data.frame(
    x,
    y
  )

  type <- match.arg(type)
  if (type == "ma") {
    if (log) {
      df[] <- lapply(df, log10)
    }
    fun <- param_plot_ma
    mu_df <- data.frame(
      colMeans(chain1@parameters[["mu"]]),
      colMeans(chain2@parameters[["mu"]])
    )
    mus <- rowMeans(mu_df)
  } else {
    fun <- param_plot_std
    mus <- NULL
  }

  colnames(df) <- c(xname, yname)
  g <- fun(
    df,
    xname = xname,
    yname = yname,
    title = paste(
      if (type == "ma") "MA plot" else "Scatter plot",
      "of HPD interval\nof", param, "parameter in",
      yname, "vs", xname, "MCMC"
    ),
    log = log && type == "std",
    mus = mus,
    ...
  )

  # ggMarginal(g)
  g
}

param_plot_ma <- function(
    df,
    xname,
    yname,
    mus,
    ...
  ) {

  x <- df[[xname]]
  y <- df[[yname]]

  df <- data.frame(
    M = y - x,
    mu = mus
  )
  m <- max(abs(df$M))
  ylim <- c(-m, m)

  g <- param_plot(df, "mu", "M", ...) +
    # geom_hline(yintercept = 0, lty = "twodash", col = "grey70", na.rm = TRUE) +
    ylim(ylim) +
    scale_x_log10()
  g
}

param_plot <- function(
    df,
    xname,
    yname,
    title,
    log = FALSE,
    bins = 100,
    ...
  ) {

  g <- ggplot(
    df,
    mapping = aes_string(
      x = paste_aes(xname),
      y = paste_aes(yname)
    )
  ) +
    # geom_point(
    #   alpha = 0.01, na.rm = TRUE) +
    # stat_density2d(
    #     aes(
    #       fill = sqrt(..density..)),
    #       geom = "tile",
    #       contour = FALSE,
    #       n = 256
    #   )  +
    scale_fill_viridis(direction = 1, name = "Density") +
    geom_hex(aes(fill = ..density..), bins = bins, na.rm = TRUE) +
    guides(fill = "none") +
    xlab(xname) +
    ylab(yname) +
    ggtitle(title)
  if (log) {
    g <- g +
      scale_x_log10() +
      scale_y_log10()
  }
  g
}

extract_hpd_interval <- function(chain, param) {
  summ <- BASiCS::Summary(chain)
  abs(summ@parameters[[param]][, 3] - summ@parameters[[param]][, 2])
}

do_hpd_plots <- function(i, maxdf, references) {
  cat(i, "/", nrow(maxdf), "\n")
  maxdf <- maxdf[i, ]
  d <- maxdf[["data"]]
  b <- maxdf[["by"]]
  rc <- references[[which(references$data == d)[[1]], "chain"]]
  c <- readRDS(maxdf[["file"]][[1]])
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
    )
  }
  c <- BASiCS:::.offset_correct(rc, c)
  l2 <- if (b == "advi") "ADVI" else "D & C"
  l3 <- sprintf(
    "log2(%s HPD interval width / Reference HPD interval width)",
    l2
  )
  l <- list(
    plot_hpd_interval(
      rc,
      "Reference",
      c,
      l2,
      "mu",
      type = "ma",
      log = FALSE,
      bins = 50
    ) + 
      labs(
        title = NULL,
        x = bquote(mu[i]),
        y = l3
      ) +
      geom_hline(aes(yintercept = 0)),
    plot_hpd_interval(
      rc,
      "Reference",
      c,
      l2,
      "epsilon",
      type = "ma",
      log = FALSE,
      bins = 50
    ) +
      labs(
        title = NULL,
        x = bquote(mu[i]),
        y = l3
      ) +
      geom_hline(aes(yintercept = 0))
  )
  ggsave(l[[1]],
    file = here(sprintf("figs/hpd/mu_%s_%s.pdf", b, d)),
    width = 7, height = 5
  )
  ggsave(l[[2]],
    file = here(sprintf("figs/hpd/epsilon_%s_%s.pdf", b, d)),
    width = 7, height = 5
  )
  l
}

do_ess_plot <- function(i, maxdf, references) {
  cat(i, "/", nrow(maxdf), "\n")
  maxdf <- maxdf[i, ]
  d <- maxdf[["data"]]
  b <- maxdf[["by"]]
  rc <- references[[which(references$data == d)[[1]], "chain"]]
  c <- readRDS(maxdf[["file"]][[1]])
  if (is.list(c)) {  
    suppressMessages(
      c <- BASiCS:::.combine_subposteriors(
        c,
        GeneOrder = rownames(rc),
        CellOrder = colnames(rc),
        SubsetBy = "gene"
      )
    )
  }
  c <- BASiCS:::.offset_correct(rc, c)

  ess <- effectiveSize(as.mcmc(c@parameters[["epsilon"]]))

  df <- data.frame(
    x = colMedians(rc@parameters[["mu"]]),
    y = colMedians(rc@parameters[["epsilon"]] - c@parameters[["epsilon"]]),
    color = ess
  )
  df <- df[order(df$color, decreasing = TRUE), ]
  ggplot(df, aes(x = x, y = y, color = color)) + 
    geom_point(alpha = 0.75) +
    scale_x_log10() +
    labs(x = bquote(mu[i]), y = bquote(epsilon[i]^Ref - epsilon[i]^Scalable)) +
    scale_color_viridis(
      name = bquote('Effective sample size'~epsilon[i]),
      trans = "log10"
    )
  ggsave(
    file = here(sprintf("figs/ess/%s_%s.pdf", b, d)),
    width = 5, height = 4
  )
}

paste_aes <- function(...) {
  paste0("`", ..., "`")
}
