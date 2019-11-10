suppressPackageStartupMessages({
  library("rstan")
  library("here")
  library("devtools")
  load_all(here("../BASiCS"))
  load_all(here())
})
source(here("data-raw/stan/functions.R"))

d <- zd()

spikes <- d[isSpike(d), ]
counts <- d[!isSpike(d), ]
spikes <- assay(spikes)
counts <- assay(counts)

L <- 12


start <- BASiCS:::HiddenBASiCS_MCMC_Start(d, 
  eta = 5, 
  m = rep(0, L), 
  V = diag(L), 
  a.sigma2 = 2, 
  b.sigma2 = 2, 
  WithSpikes = TRUE)

ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)
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
  mbeta = rep(0, L),
  ml = ml,
  l = 12,
  vbeta = diag(L),
  rbf_variance = 1.2,
  eta = 5,
  astwo = 2,
  bstwo = 2
)



t <- system.time({
  set.seed(42)
  b <- BASiCS_MCMC(
    d,
    ml  = ml,
    N = 20000, 
    Thin = 10, 
    Burn = 10000, 
    Regression = TRUE, 
    WithSpikes = TRUE)
})
saveRDS(t, file = here("outputs/basics_time.rds"))
saveRDS(b, file = here("outputs/basics_zeisel.rds"))




g1 <- BASiCS_showFit(b_reg)
g2 <- BASiCS_showFit(b_reg, ml = sdata$ml)



t <- system.time({
  set.seed(42)
  v <- vb(stanmodels$basics_regression, iter = 1000000, data = sdata)    
})
saveRDS(t, file = here("outputs/vb_time.rds"))
saveRDS(v, file = here("outputs/vb_zeisel.rds"))

t <- system.time({
  set.seed(42)
  m <- sampling(
    stanmodels$basics_regression, 
    data = sdata)  
})

saveRDS(t, file = here("outputs/hmc_time.rds"))
saveRDS(m, file = here("outputs/hmc_zeisel.rds"))
rm(m)

