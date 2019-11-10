stan_plot <- function(chain, param = "mu", which = 1:10) {
  chain <- Summary(chain)
  dat <- chain@parameters[[param]][which, , drop = FALSE]
  dat <- as.data.frame(dat)
  dat$ID <- factor(rownames(dat), levels=rev(rownames(dat)))
  ggplot(dat, aes(x=ID, y=median, ymin=lower, ymax=upper)) + 
    geom_pointrange() +
    coord_flip()
}


hpd_plot <- function(c1, c2, param = "mu") {
  p1 <- extract_param(c1, param = param)
  p2 <- extract_param(c2, param = param)
  if (param != "mu") {
    m1 <- extract_param(c1, param = "mu")
    m2 <- extract_param(c2, param = "mu")
  } else {
    m1 <- p1
    m2 <- p2
  }
  p1$width <- hpd_width(p1)
  p2$width <- hpd_width(p2)
  df <- p1
  df$width_diff <- log2(abs(p1$width / p2$width))
  df$mean_mu <- rowMeans(data.frame(m1$median, m2$median))

  ggplot(df, aes(x = mean_mu, y = width_diff)) + 
    labs(x = "Mean expression level", y = "log2(Ratio of HPD interval width)") +
    # geom_hline(yintercept = c(1, -1), colour = "grey60", lty = "dashed") +
    geom_hline(yintercept = 0) +
    geom_hex(bins = 100) +
    ylim(-max(abs(df$width_diff)), max(abs(df$width_diff))) +
    scale_x_log10() +
    scale_fill_viridis(name = NULL)
}

hpd_width <- function(df) {
  abs(df$upper - df$lower)
}

ma_plot <- function(c1, c2, param = "mu", log = TRUE) {
  p1 <- extract_param(c1, param = param)
  p2 <- extract_param(c2, param = param)
  if (param != "mu") {
    m1 <- extract_param(c1, param = "mu")
    m2 <- extract_param(c2, param = "mu")
  } else {
    m1 <- p1
    m2 <- p2
  }
  df <- p1
  df$a <- rowMeans(data.frame(m1$median, m2$median))
  df$m <- p1$median / p2$median
  if (log) df$m <- log2(df$m)

  ggplot(df, aes(x = a, y = m)) + 
    labs(x = "Mean expression level", y = "log2(ratio)") +
    geom_hline(yintercept = c(-log2(1.5), log2(1.5)), colour = "grey60", lty = "dashed") +
    geom_hline(yintercept = 0) +
    geom_hex(bins = 100) +
    ylim(-max(abs(df$m)), max(abs(df$m))) +
    scale_x_log10() +
    scale_fill_viridis(name = NULL)
}
