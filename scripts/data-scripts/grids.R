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

advi_grid <- expand.grid(datasets, seeds + 1)

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

divide_and_conquer_grid <- expand.grid(datasets, chains, seeds, subsets, stringsAsFactors = FALSE)
colnames(divide_and_conquer_grid) <- c("dataset", "chains", "seed", "subset")

divide_and_conquer_grid <- merge(divide_and_conquer_grid, dims)

ind_gene <- divide_and_conquer_grid$nGene / divide_and_conquer_grid$chains > 150
divide_and_conquer_grid <- divide_and_conquer_grid[ind_gene, ]

ind_one <- divide_and_conquer_grid$chain > 1 | 
  (divide_and_conquer_grid$chain == 1 & divide_and_conquer_grid$seed == 7)
divide_and_conquer_grid <- divide_and_conquer_grid[ind_one, ]


write.table(
  divide_and_conquer_grid[, c("dataset", "chains", "seed", "subset")],
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/divide_and_conquer_grid.txt"
)


downsampling_grid <- expand.grid("zeisel", seq(1, 0.2, length.out = 5), seeds)
ind_one <- downsampling_grid$Var2 == 1 & downsampling_grid$Var3 == 7 |
  downsampling_grid$Var2 != 1
downsampling_grid <- downsampling_grid[ind_one, ]
write.table(
  downsampling_grid,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/downsampling_grid.txt"
)
