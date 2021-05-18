library("BASiCS")

n <- 200
q <- 5000
p <- 0.2
fc <- 2

Mu1 <- rlnorm(q)
Delta1 <- rlnorm(q, (log(Mu1) * -1) + 1)
S1 <- rgamma(n, 1, 1)
Theta1 <- rgamma(2, 1, 1)
BatchInfo1 <- rbinom(n, 1, 0.5)

sce1 <- BASiCS_Sim(
  Mu = Mu1,
  Delta = Delta1,
  S = S1,
  Theta = Theta1,
  BatchInfo = BatchInfo1
)

IndDE <- rbinom(q, 1, p)
IndDD <- rbinom(q, 1, p)
Mu2 <- exp(log(Mu1) + (sample(c(-fc, fc), q, replace = TRUE) * IndDE))
Delta2 <- exp(log(Delta1) + (sample(c(-fc, fc), q, replace = TRUE) * IndDD))
S2 <- rgamma(n, 1, 1)
Theta2 <- rgamma(2, 1, 1)
BatchInfo2 <- rbinom(n, 1, 0.5)

sce2 <- BASiCS_Sim(
  Mu = Mu2,
  Delta = Delta2,
  S = S2,
  Theta = Theta2,
  BatchInfo = BatchInfo2
)

saveRDS(
  list(
    sces=list(sce1, sce2),
    IndDE = IndDE,
    IndDD = IndDD,
    Mus = list(Mu1, Mu2),
    Deltas = list(Delta1, Delta2),
    S = list(S1, S2),
    Thetas = list(Theta1, Theta2),
    BatchInfo = list(BatchInfo1, BatchInfo1)
  ),
  "data/simulation.rds"
)
