library("SingleCellExperiment")
library("here")

source(here("scripts/chain-scripts/benchmark_code.R"))

seeds <- benchmarkSeeds()

datasets <- c(
  "tung",
  "buettner",
  "pbmc",
  # "splatter", 
  "zeisel"
  # , "williams"
)

advi_grid <- expand.grid(datasets, seeds)

write.table(
  advi_grid,
  file = "data/advi_grid.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)


chains <- c(
  1, 
  2, 
  4, 
  8, 
  16, 
  32,
  64, 
  128,
  256
)


subsets <- c(
  # "cell", 
  "gene"
)

dims <- sapply(datasets, function(x) {
  dim(readRDS(paste0("data/", x, ".rds")))
})
rownames(dims) <- c("nGene", "nCell")
dims <- as.data.frame(t(dims))
dims$dataset <- rownames(dims)

grid <- expand.grid(datasets, chains, seeds, subsets, stringsAsFactors = FALSE)
colnames(grid) <- c("dataset", "chains", "seed", "subset")

grid <- merge(grid, dims)

grid <- grid[grid$nGene / grid$chains > 150, ]

grid <- grid[grid$chain > 1 | (grid$chain == 1 & grid$seed == 7), ]


write.table(
  grid[, c("dataset", "chains", "seed", "subset")],
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/divide_and_conquer_grid.txt"
)
