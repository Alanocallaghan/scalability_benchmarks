suppressPackageStartupMessages({
  library("rstan")
  library("devtools")
  load_all("../BASiCS")
  library("here")
  load_all()
})
source(here("data-raw/stan/functions.R"))
options(mc.cores = 4)
theme_set(theme_bw())

set.seed(42)
d <- makeExampleBASiCS_Data(WithSpikes = TRUE)

spikes <- d[isSpike(d), ]
counts <- d[!isSpike(d), ]
spikes <- assay(spikes)
counts <- assay(counts)

l <- 12

start <- BASiCS:::HiddenBASiCS_MCMC_Start(d,
  eta = 5,
  m = rep(0, l),
  V = diag(l),
  a.sigma2 = 2,
  b.sigma2 = 2,
  WithSpikes = TRUE
)
ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, l)

sdata <- list(
    g = nrow(counts), 
    c = ncol(counts),
    sg = nrow(spikes),
    counts = counts, 
    spikes = spikes,
    spike_levels = metadata(d)$SpikeInput,
    as = 1,
    bs = 1,
    atheta = 1,
    btheta = 1,
    smu = sqrt(0.5),
    sdelta = sqrt(0.5),
    aphi = rep(1, ncol(counts)),
    mbeta = rep(0, l),
    l = l,
    ml = ml,
    vbeta = diag(l),
    rbf_variance = 1.2,
    eta = 5,
    astwo = 2,
    bstwo = 2
)

set.seed(42)
b <- BASiCS_MCMC(
  d,
  N = 20000, 
  Thin = 10, 
  Burn = 10000, 
  Regression = FALSE, 
  PrintProgress = FALSE,
  WithSpikes = TRUE
)

set.seed(42)
b_reg <- BASiCS_MCMC(
  d,
  N = 20000, 
  Thin = 10, 
  Burn = 10000, 
  ml = sdata$ml,
  Regression = TRUE, 
  PrintProgress = FALSE,
  WithSpikes = TRUE
)


set.seed(42)
m <- sampling(stanmodels$basics, data = sdata)

set.seed(42)
v <- vb(stanmodels$basics, iter = 100000, data = sdata)    


set.seed(42)
m_reg <- sampling(stanmodels$basics_regression, data = sdata)

set.seed(42)
system.time(
c <- capture.output(v_reg <- vb(stanmodels$basics_regression, iter = 100000, data = sdata))
)


save.image("basics_stan_example.RData")


g1 <- comp_plot(b_reg, "AMWG", m_reg, "HMC", v_reg, "ADVI", "beta")


g2 <- comp_plot(b_reg, "BASiCS", m_reg, "HMC", v_reg, "ADVI", "mu")
g3 <- comp_plot(b_reg, "BASiCS", m_reg, "HMC", v_reg, "ADVI", "delta")
g4 <- comp_plot(b_reg, "BASiCS", m_reg, "HMC", v_reg, "ADVI", "epsilon")


ggsave(
  g1 + theme(legend.position = "bottom"),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/stan_mu.pdf",
  width = 4,
  height = 4
)



ggsave(
  g1 + theme(legend.position = "bottom"),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/beta_hmc_advi.pdf",
  width = 6, 
  height = 4
)
ggsave(
  g2 + theme(legend.position = "bottom"),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/mu_stan.pdf",
  width = 8, 
  height = 6
)
ggsave(
  g3 + theme(legend.position="bottom"),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/delta_stan.pdf",
  width = 8, 
  height = 6
)
ggsave(
  g4 + theme(legend.position = "bottom"),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/epsilon_stan.pdf",
  width = 8, 
  height = 6
)




# vs <- extract(v_reg)
# vc <- new("BASiCS_Chain", 
#   parameters = list(
#     mu = vs$mu, 
#     delta = vs$delta, 
#     nu = vs$nu,
#     s = vs$s,
#     epsilon = vs$epsilon,
#     beta = vs$beta,
#     phi = vs$phi,
#     theta = vs$theta
#   )
# )

# mc <- new("BASiCS_Chain", 
#   parameters = list(
#     mu = ms$mu, 
#     delta = ms$delta, 
#     nu = ms$nu,
#     s = ms$s,
#     phi = ms$phi,
#     theta = ms$theta
#   )
# )
