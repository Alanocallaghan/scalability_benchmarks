library("SingleCellExperiment")
library("ggplot2")

theme_set(theme_bw())
# lapply(
#   c(
#     "tung",
#     "buettner",
#     "pbmc",
#     # "splatter",
#     "zeisel"
#     # , "williams"
#   ),
#   function(d) table(readRDS(paste0("data/", d, ".rds"))@colData$BatchInfo)
# )

datasets <- c(
  "tung",
  "buettner",
  # "pbmc",
  "chen",
  "zeisel"
)
datas <- lapply(datasets, function(d) readRDS(paste0("rdata/", d, ".rds")))

names(datas) <- c(
  "Tung",
  "Buettner",
  # "10x PBMC",
  "Chen",
  "Zeisel"
)

l <- lapply(names(datas),
  function(n) {
    x <- datas[[n]]
    data.frame(
      name = n,
      lib_sizes = colSums(counts(x)),
      num_feats = colSums(counts(x) != 0)
    )
  }
)
data_df <- do.call(rbind, l)

g <- ggplot(data_df, aes(x = lib_sizes, color = name, fill = name)) +
  labs(x = "Library size") +
  scale_color_brewer(
    palette = "Set2",
    name = "Dataset",
    aesthetics = c("color", "fill")
  ) +
  geom_density(alpha = 0.2) +
  scale_x_log10()

ggsave("figs/libsize_density.pdf", width = 5, height = 4)

g <- ggplot(data_df, aes(x = num_feats, color = name, fill = name)) +
  labs(x = "Number of expressed features") +
  scale_color_brewer(
    palette = "Set2",
    name = "Dataset",
    aesthetics = c("color", "fill")
  ) +
  geom_density(alpha = 0.2) +
  scale_x_log10()
ggsave("figs/complexity_density.pdf", width = 5, height = 4)

l <- lapply(names(datas),
  function(n) {
    x <- datas[[n]]
    data.frame(
      name = n,
      mean_expression = rowMeans(counts(x)),
      dropout = rowMeans(counts(x) == 0)
    )
  }
)
data_df <- do.call(rbind, l)

g <- ggplot(data_df, aes(x = mean_expression, color = name, fill = name)) +
  labs(x = "Mean expression") +
  scale_color_brewer(
    palette = "Set2",
    name = "Dataset",
    aesthetics = c("color", "fill")
  ) +
  geom_density(alpha = 0.2) +
  scale_x_log10()
ggsave("figs/expression_density.pdf", width = 5, height = 4)

g <- ggplot(data_df, aes(x = dropout, color = name, fill = name)) +
  labs(x = "Proportion of zeros") +
  scale_color_brewer(
    palette = "Set2",
    name = "Dataset",
    aesthetics = c("color", "fill")
  ) +
  geom_density(alpha = 0.2)
ggsave("figs/dropout_density.pdf", width = 5, height = 4)
