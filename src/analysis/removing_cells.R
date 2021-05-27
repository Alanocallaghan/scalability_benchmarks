library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")

theme_set(theme_bw())
source(here("src/analysis/functions.R"))

rm_files <- list.files("outputs/removing/divide", full.names = TRUE)
rmt <- file2triplets(rm_files)
rmt <- rmt[as.logical(sapply(rmt, length))]
rm_df <- read_triplets(rmt, combine = TRUE)

rm_ref_files <- list.files("outputs/removing/reference", full.names = TRUE)
ref_df_rm <- read_triplets(file2triplets(rm_ref_files), combine = TRUE)
ref_df_rm$chain <- lapply(ref_df_rm$file, readRDS)

rm_df <- do_de(rm_df, ref_df = ref_df_rm, match_column = "proportion_retained")

mdf_rm <- reshape2::melt(
  rm_df,
  measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_rm$variable <- gsub("pDiffExp", "mu", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffDisp", "delta", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffResDisp", "epsilon", mdf_rm$variable)
mdf_rm$variable <- factor(mdf_rm$variable, levels = c("mu", "delta", "epsilon"))
mdf_rm$cells_retained <- mdf_rm$proportion_retained * mdf_rm$nCells

mdf_rm$proportion_retained <- factor(
  paste(mdf_rm$proportion_retained * 100, "%"),
  levels = paste(sort(unique(rm_df$proportion_retained), decreasing = TRUE) * 100, "%")
)

ggplot(mdf_rm, aes(x = cells_retained, y = value, color = variable)) +
  geom_quasirandom(dodge.width = 60, size = 0.4) +
  scale_color_brewer(name = "Parameter", palette = "Set1") +
  scale_x_reverse("Number of cells") +
  scale_y_continuous(label = scales::percent) +
  labs(x = "nCells", y = "Portion of genes perturbed")

ggsave("figs/removing_cells.pdf", width = 6, height = 4)
