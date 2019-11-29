options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())
source(here("scripts/analysis-scripts/functions.R"))

ds_files <- list.files("outputs/downsampling/divide/", full.names = TRUE)
dst <- file2triplets(ds_files)
dst <- dst[as.logical(sapply(dst, length))]
ds_df <- read_triplets(dst, combine = TRUE)

ref_files <- list.files("outputs/downsampling/reference", full.names = TRUE)
ref_df_ds <- read_triplets(file2triplets(ref_files), combine = TRUE)
ref_df_ds$chain <- lapply(ref_df_ds$file, readRDS)
ds_df <- do_de(ds_df, ref_df_ds, "downsample_rate")

mdf_ds <- reshape2::melt(ds_df,
  measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_ds$variable <- gsub("pDiffExp", "mu", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffDisp", "delta", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffResDisp", "epsilon", mdf_ds$variable)
mdf_ds$variable <- factor(mdf_ds$variable, levels = c("mu", "delta", "epsilon"))
mdf_ds$downsample_rate <- factor(
  paste(mdf_ds$downsample_rate * 100, "%"),
  levels = paste(sort(unique(ds_df$downsample_rate), decreasing = TRUE) * 100, "%")
)

ggplot(mdf_ds, aes(x = downsample_rate, y = value, color = variable)) +
  geom_quasirandom(dodge.width = 0.5, size = 0.25) +
  scale_color_brewer(name = "Parameter", palette = "Set1") +
  # scale_x_continuous(label = scales::percent) +
  scale_y_continuous(label = scales::percent) +
  labs(x = "Proportion of original counts", y = "Portion of genes perturbed")

ggsave("figs/downsampling.pdf", width = 6, height = 4)
