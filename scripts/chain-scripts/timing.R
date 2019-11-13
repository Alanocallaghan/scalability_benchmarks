library("Scalability")
library("BASiCS")

args <- commandArgs(trailingOnly = TRUE)
dir.create("outputs/time/", recursive = TRUE, showWarnings = FALSE)
data <- args[[1]]

time_mcmc <- function(n, times = 1) {
  replicate(times, {
    subsets <- Scalability:::generateSubsets(
    	readRDS(paste0("data/", data, ".rds")),
    	NSubsets = n,
    	SubsetBy = "gene",
    	WithSpikes = TRUE
    )
    system.time(
      suppressMessages(
      	capture.output(
	        BASiCS_MCMC(
	          subsets[[1]],
	          N = 20000,
	          Thin = 10,
	          Burn = 10000,
	          WithSpikes = TRUE,
	          Regression = TRUE,
	          PrintProgress = FALSE
	        )
    		)
      )
    )[["elapsed"]]
  })
}

for (n in c(2, 4, 8, 16, 32)) {
  time <- time_mcmc(n, times = 6)
  saveRDS(time, paste0("outputs/time/", data, "_", n, ".rds"))
}
