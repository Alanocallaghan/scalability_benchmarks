## Runs BASiCS in a divide and conquer setting
divide_and_conquer_benchmark <- function(
        Data,
        N = 20000,
        Thin = 10,
        Burn = 10000,
        DataName,
        NSubsets = 1,
        SubsetBy = c("gene", "cell"),
        Regression = TRUE,
        WithSpikes = length(SingleCellExperiment::altExpNames(Data)) > 0,
        Seed = 42,
        ...) {
    set.seed(Seed)
    SubsetBy <- match.arg(SubsetBy)
    t <- system.time(
        chain <- BASiCS_MCMC(
            Data,
            NSubsets = NSubsets,
            SubsetBy = SubsetBy,
            N = N,
            Thin = Thin,
            Burn = Burn,
            WithSpikes = WithSpikes,
            PrintProgress = FALSE,
            Regression = Regression,
            ...
        )
    )
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
