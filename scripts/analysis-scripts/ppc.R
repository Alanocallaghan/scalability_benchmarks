datas <- sapply(
  c("tung", "zeisel", "buettner", "pbmc"),
  function(dataset) {
    readRDS(paste0("data/", dataset, ".rds"))
  },
  simplify = FALSE
)

## distribution of mean expression
get_means <- function(datas) {
  sapply(datas, function(x) rowMeans(counts(x)), simplify = FALSE)
}

## distribution of proportion of zeros
get_dropout <- function(datas) {
  sapply(datas, function(x) rowMeans(counts(x) != 0), simplify = FALSE)
}

## distribution of cell size
get_size <- function(datas) {
  sapply(datas, function(x) colSums(counts(x)), simplify = FALSE)
}

## distribution of cell complexity
get_complexity <- function(datas) {
  sapply(datas, function(x) colSums(counts(x) != 0), simplify = FALSE)
}


base_means <- get_means(datas)
base_dropout <- get_dropout(datas)
base_size <- get_size(datas)
base_complexity <- get_complexity(datas)


file <- list.files("outputs/divide_and_conquer", full.names = TRUE)[20]
chains <- readRDS(paste0(file, "/chains.rds"))
config <- readRDS(paste0(file, "/config.rds"))
if (is.list(chains)) {
  chains <- combine_subposteriors(
    chains,
    subset_by = config$by
  )
}

ppc_datas <- lapply(1:50, function(i) {
  BASiCS_Draw(chains, N = i)
})

ppc_means <- get_means(ppc_datas)
ppc_dropout <- get_dropout(ppc_datas)
ppc_size <- get_size(ppc_datas)
ppc_complexity <- get_complexity(ppc_datas)


density_plot <- function(base_param, ppc_param) {
  df <- do.call(cbind, ppc_param)
  df <- as.data.frame(df)
  df <- cbind(df, base = base_param[[config$data]])
  mdf <- reshape2::melt(as.matrix(df))
  ggplot(mdf, aes(x = value, group = Var2, color = Var2 == "base")) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Base", "Posterior sample"),
      guide = FALSE
    ) +
    geom_density() +
    scale_x_log10()  
}

density_plot(base_means, ppc_means)
density_plot(base_dropout, ppc_dropout)
density_plot(base_size, ppc_size)
density_plot(base_complexity, ppc_complexity)


sapply(1:50, 
  function(i) {
    ks.test(ppc_means[[i]], base_means[[config$data]])$p.value
  }
)

