suppressPackageStartupMessages({
  library("ggplot2")
  library("rstan")
  library("BASiCS")
  library("viridis")
})
theme_set(theme_bw())
source("functions.R")


b <- readRDS("basics_zeisel.rds")
v <- readRDS("vb_zeisel.rds")
# m <- readRDS("hmc_zeisel.rds")



g <- ma_plot(b, v, "mu")
ggsave(g, file = "basics_advi_ma_mu.pdf", width = 6, height = 4)
g <- ma_plot(b, v, "delta")
ggsave(g, file = "basics_advi_ma_delta.pdf", width = 6, height = 4)
ma_plot(b, v, "epsilon")



g <- hpd_plot(b, v, "mu")
ggsave(g, file = "basics_advi_hpd_mu.pdf", width = 6, height = 4)
g <- hpd_plot(b, v, "delta")
ggsave(g, file = "basics_advi_hpd_delta.pdf", width = 6, height = 4)
hpd_plot(b, v, "epsilon")
system("cp *.pdf /home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/")



vs <- extract(v)

params <- c("mu", "delta", "nu", "s", "phi", "theta")
vs[params] <- lapply(params, function(name) {
  if (length(dim(vs[[name]])) > 1) {
    colnames(vs[[name]]) <- colnames(b@parameters[[name]])
  }
  vs[[name]]
})

vc <- new("BASiCS_Chain", 
  parameters = list(
    mu = vs$mu, 
    delta = vs$delta, 
    nu = vs$nu,
    s = vs$s,
    phi = vs$phi,
    theta = vs$theta
  )
)
de <- BASiCS_TestDE(
  b, 
  vc, 
  GroupLabel1 = "BASiCS", 
  GroupLabel2 = "ADVI", 
  Plot = FALSE, 
  PlotOffset = FALSE)

