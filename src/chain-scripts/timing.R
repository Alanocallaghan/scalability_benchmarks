suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
options(stringsAsFactors=FALSE)
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-f", "--fraction", type = "double")
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
            WithSpikes = data != "pbmc",
            Regression = TRUE,
            PrintProgress = FALSE
          )
        )
      )
    )[["elapsed"]]
  })
}

for (n in c(2, 4, 8, 16, 32)) {
  cat(n, "chains\n")
  time <- time_mcmc(n, times = 6)
  saveRDS(time, paste0("outputs/time/", data, "_", n, ".rds"))
}
