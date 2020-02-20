#' Random seed used for benchmarking
benchmarkSeeds <- function() {
  c(7, 
    14, 
    21, 
    28, 
    35, 
    42
  )
}


## Runs BASiCS in a divide and conquer setting
divide_and_conquer_benchmark <- function(
    Data,
    N = 20000,
    Thin = 10,
    Burn = 10000,
    DataName,
    NSubsets = 1,
    SubsetBy = c("gene", "cell"),
    Cores = min(NSubsets, detectCores(all.tests = FALSE, logical = TRUE)),
    Regression = TRUE,
    WithSpikes = length(SingleCellExperiment::altExpNames(Data)) > 0,
    Seed = 42,
    ...) {

  SubsetBy <- match.arg(SubsetBy)
  if (NSubsets == 1) {
    t <- system.time(
      chain <- BASiCS_MCMC(
        Data,
        N = N,
        Thin = Thin,
        Burn = Burn,
        WithSpikes = WithSpikes,
        PrintProgress = FALSE,
        Regression = Regression,
        ...
      )
    )
  } else {
      try({
        t <- system.time(
          chain <- Scalability:::multi_MCMC(
            Data,
            NSubsets = NSubsets,
            mc.cores = Cores,
            SubsetBy = SubsetBy,
            N = N,
            Thin = Thin,
            Burn = Burn,
            WithSpikes = WithSpikes,
            PrintProgress = FALSE,
            Regression = Regression,
            Seed = Seed,
            ...
          )
      )
    })
  }
  list(
    time = t,
    chain = chain, 
    config = list(
      chains = NSubsets,
      by = SubsetBy,
      data = DataName,
      seed = Seed
    )
  )
}
