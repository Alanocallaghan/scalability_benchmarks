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
cat("Doing", data, "\n")
time_mcmc <- function(n, times = 1) {
  replicate(times, {
    subsets <- BASiCS:::.generateSubsets(
      readRDS(paste0("data/", data, ".rds")),
      NSubsets = n,
      SubsetBy = "gene",
      WithSpikes = data != "pbmc"
    )
    system.time(
      suppressMessages(
        capture.output(
          BASiCS_MCMC(
            subsets[[1]],
            N = args[["iterations"]],
            Thin = max((args[["iterations"]] / 2) / 1000, 2),
            Burn = max(args[["iterations"]] / 2, 4),
            WithSpikes = "spike-ins" %in% altExpNames(data),
            Regression = TRUE,
            PrintProgress = FALSE
          )
        )
      )
    )[["elapsed"]]
  })
}
n <- args[["nsubsets"]]
time <- time_mcmc(n, times = 6)
saveRDS(time, args[["output"]])
