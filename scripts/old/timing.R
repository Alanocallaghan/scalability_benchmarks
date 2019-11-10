devtools::load_all()

d <- zd()


time_mcmc <- function(n, times = 1) {
  replicate(times, {
    ds <- generateSubsets(d, NSubsets = n, SubsetBy = "gene", WithSpikes = TRUE)
#    setMKLthreads(1)
    times <- mclapply(seq_along(ds),
      function(i) {
        cat(i, "\n")
        system.time(
          suppressMessages(capture.output({
            BASiCS_MCMC(
              ds[[i]],
              N = 20000,
              Thin = 10,
              Burn = 10000,
              WithSpikes = TRUE,
              Regression = TRUE,
              PrintProgress = FALSE
            )
          }
        ))
        )[["elapsed"]]
      },
      mc.cores = 16
    )
    unlist(times)
  })
}

t32 <- time_mcmc(32, times = 2)
t64 <- time_mcmc(64, times = 1)
t128 <- time_mcmc(128, times = 1)

saveRDS(t32, "outputs/t32.rds")
saveRDS(t64, "outputs/t128.rds")
saveRDS(t128, "outputs/t128.rds")
