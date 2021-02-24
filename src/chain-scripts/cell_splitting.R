library("BASiCS")

theme_set(theme_bw())

fit <- multi_MCMC(
  readRDS("data/zeisel.rds"),
  SubsetBy = "cell",
  NSubsets = 32,
  N = 20000,
  Thin = 10,
  Burn = 10000
)
dir.create("outputs/cell_splitting")
saveRDS(fit, "outputs/cell_splitting/zeisel.rds")
ref_file <- (df %>% filter(chains == 1, data == "zeisel") %>% pull(file))[[1]]

ref <- readRDS(ref_file)

fitc <- BASiCS:::.combine_subposteriors(
  fit,
  GeneOrder = colnames(ref@parameters[["mu"]]),
  CellOrder = colnames(ref@parameters[["nu"]]),
  SubsetBy = "cell"
)

d <- BASiCS_TestDE(
  fitc,
  ref,
  GroupLabel1 = "D&C",
  GroupLabel2 = "Reference"
)

g <- BASiCS_PlotDE(
  d@Results[[1]],
  Plots = c("MAPlot")
)

ggsave(g, file = "figs/cell_partitions.pdf", width = 6, height = 4)
