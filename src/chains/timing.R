suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-n", "--nsubsets", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()


dir.create("outputs/time/", recursive = TRUE, showWarnings = FALSE)
data <- args[["data"]]
# cat("Doing", data, "\n")
sce <- readRDS(paste0("data/", data, ".rds"))
spikes <- "spike-ins" %in% altExpNames(sce)
time_mcmc <- function(n, times = 1) {
  if (n == 1) {
    replicate(times, {
      system.time(
        suppressMessages(
          capture.output(
            BASiCS_MCMC(
              sce,
              N = args[["iterations"]],
              Thin = max((args[["iterations"]] / 2) / 1000, 2),
              Burn = max(args[["iterations"]] / 2, 4),
              WithSpikes = spikes,
              Regression = TRUE,
              PrintProgress = FALSE
            )
          )
        )
      )[["elapsed"]]
    })
  } else {
    replicate(times, {
      subsets <- BASiCS:::.generateSubsets(
        sce,
        NSubsets = n,
        SubsetBy = "gene",
        WithSpikes = spikes
      )
      system.time(
        suppressMessages(
          capture.output(
            BASiCS_MCMC(
              subsets[[1]],
              N = args[["iterations"]],
              Thin = max((args[["iterations"]] / 2) / 1000, 2),
              Burn = max(args[["iterations"]] / 2, 4),
              WithSpikes = spikes,
              Regression = TRUE,
              PrintProgress = FALSE
            )
          )
        )
      )[["elapsed"]]
    })
  }
}
n <- args[["nsubsets"]]
time <- time_mcmc(n, times = 6)
saveRDS(time, args[["output"]])
