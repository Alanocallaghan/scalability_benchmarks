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

advi_grid <- expand.grid(datasets, seeds - 1)

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


removing_grid <- expand.grid("zeisel", seq(1, 0.1, length.out = 10), seeds)
removing_grid_ref <- expand.grid("zeisel", seq(1, 0.1, length.out = 10))

write.table(
  removing_grid,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/removing_grid.txt"
)

write.table(
  removing_grid_ref,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/removing_grid_ref.txt"
)


downsampling_grid <- expand.grid("buettner", seq(1, 0.1, length.out = 10), seeds)
downsampling_ref_grid <- expand.grid("buettner", seq(1, 0.1, length.out = 10))

write.table(
  downsampling_grid,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/downsampling_grid.txt"
)

write.table(
  downsampling_ref_grid,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/downsampling_grid_ref.txt"
)






write.table(
  data.frame(datasets),
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/datasets.txt"
)


write.table(
  data.frame(c(
    "tung",
    # "buettner",
    # "pbmc",
    # "splatter", 
    "zeisel"
    # , "williams"
    )
  ),
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  file = "data/datasets_batch.txt"
)
