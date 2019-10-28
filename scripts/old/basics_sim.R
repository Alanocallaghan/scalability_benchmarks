devtools::load_all()

d <- zeiselData()
set.seed(42)

b <- BASiCS_MCMC(
  d,
  N = 10,
  Thin = 2,
  Burn = 4,
  Regression = TRUE,
  WithSpikes = TRUE
)

ests <- lapply(b@parameters, colMedians)
s <- ests$s + rnorm(length(ests$nu) * 10, 1, 0.2)
phi <- gtools::rdirichlet(1, s)[1, ] * length(s)

sim <- BASiCS_Sim(
  Mu = ests$mu,
  Mu_spikes = rowMeans(counts(d)[isSpike(d), ]),
  Delta = ests$delta,
  Phi = phi,
  S = s,
  Theta = ests$theta
)


c <- counts(sim)

N <- 20000
Thin <- 10
Burn <- 10000

b_full <- BASiCS_MCMC(
  sim,
  N = N,
  Thin = Thin,
  Burn = Burn,
  Regression = TRUE,
  WithSpikes = TRUE
)

b_sub <- lapply(c(0.8, 0.6, 0.4, 0.2),
  function(n) {
    counts(sim) <- apply(c, 2, function(col) {
      rbinom(length(col), col, n)
    })
    BASiCS_MCMC(
      sim,
      N = N,
      Thin = Thin,
      Burn = Burn,
      Regression = TRUE,
      WithSpikes = TRUE
    )
  }
)
