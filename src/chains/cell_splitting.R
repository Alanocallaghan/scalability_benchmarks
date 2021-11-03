#!/usr/bin/env Rscript

## TODO: UNUSED

suppressPackageStartupMessages({
  library("argparse")
  library("here")
  library("BASiCS")
})
parser <- ArgumentParser()
parser$add_argument("-d", "--data")
parser$add_argument("-s", "--seed", type = "double")
parser$add_argument("-i", "--iterations", type = "double")
parser$add_argument("-o", "--output")
args <- parser$parse_args()

theme_set(theme_bw())

fit <- BASiCS_MCMC(
  readRDS("data/zeisel.rds"),
  SubsetBy = "cell",
  NSubsets = 32,
  PrintProgress = FALSE,
  N = args[["iterations"]],
  Thin = max((args[["iterations"]] / 2) / 1000, 2),
  Burn = max(args[["iterations"]] / 2, 4)
)
dir.create(args[["output"]])
saveRDS(fit, "outputs/cell_splitting/", args[["data"]], ".rds")


# ref_file <- (df %>% filter(chains == 1, data == "zeisel") %>% pull(file))[[1]]

# ref <- readRDS(ref_file)


# d <- BASiCS_TestDE(
#   fitc,
#   ref,
#   GroupLabel1 = "D&C",
#   GroupLabel2 = "Reference"
# )

# g <- BASiCS_PlotDE(
#   d@Results[[1]],
#   Plots = c("MAPlot")
# )

# ggsave(g, file = "figs/cell_partitions.pdf", width = 6, height = 4)
