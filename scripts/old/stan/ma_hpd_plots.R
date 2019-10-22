load_all()

b <- readRDS("outputs/basics_zeisel.rds")
v <- readRDS("outputs/vb_zeisel.rds")


tobasics <- function(stanfit) {
  s <- extract(stanfit)
  parameters <- list(
    mu = s$mu,
    delta = s$delta,
    epsilon = s$epsilon,
    s = s$s,
    nu = s$nu,
    theta = as.matrix(s$theta),
    phi = s$phi
  )
  gp <- c("mu", "delta", "epsilon")
  cp <- c("s", "nu", "phi")
  parameters[gp] <- lapply(gp, function(x) {
    colnames(parameters[[x]]) <- colnames(b@parameters[["mu"]])
    parameters[[x]]
  })
  parameters[cp] <- lapply(cp, function(x) {
    colnames(parameters[[x]]) <- colnames(b@parameters[["nu"]])
    parameters[[x]]
  })
  new("BASiCS_Chain", parameters = parameters)
}
bvb <- tobasics(v)
bs <- Summary(b)
bvbs <- Summary(bvb)

theme_set(theme_bw())



g <- plot_param_ma(bs, "BASiCS", bvbs, "ADVI", "mu") + 
  ggtitle(NULL) +
  labs(x = bquote(mu[i]), y = bquote(log2(mu[i]^MCMC / mu[i]^ADVI))) 
ggsave(g, 
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_ma_mu.pdf",
  width = 3, height = 3)
g <- plot_param_ma(bs, "BASiCS", bvbs, "ADVI", "delta") + 
  ggtitle(NULL) + 
  labs(x = bquote(mu[i]), y = bquote(log2(delta[i]^MCMC / delta[i]^ADVI)))
ggsave(g, 
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_ma_delta.pdf",
  width = 3, height = 3)
g <- plot_param_ma(bs, "BASiCS", bvbs, "ADVI", "epsilon", log = FALSE) + 
  ggtitle(NULL) +
  labs(x = bquote(mu[i]), y = bquote(epsilon[i]^MCMC - epsilon[i]^ADVI))
ggsave(g, 
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_ma_epsilon.pdf",
  width = 3, height = 3)




# g <- plot_hpd_interval(bs, "BASiCS", bvbs, "ADVI", "mu") + ggtitle(NULL)
# ggsave(g, 
#   file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_hpd_mu.pdf",
#   width = 4, height = 3)
# g <- plot_hpd_interval(bs, "BASiCS", bvbs, "ADVI", "delta") + ggtitle(NULL)
# ggsave(g, 
#   file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_hpd_delta.pdf",
#   width = 4, height = 3)
# g <- plot_hpd_interval(bs, "BASiCS", bvbs, "ADVI", "epsilon", log = FALSE) + ggtitle(NULL)
# ggsave(g, 
#   file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/basics_advi_hpd_epsilon.pdf",
#   width = 4, height = 3)
