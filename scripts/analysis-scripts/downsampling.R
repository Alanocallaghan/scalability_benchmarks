options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())
source(here("scripts/analysis-scripts/functions.R"))

ds_files <- list.files("/home/alan/Documents/scratchdir/downsampling/divide/", full.names = TRUE)
dst <- file2triplets(ds_files)
dst <- dst[as.logical(sapply(dst, length))]
ds_df <- read_triplets(dst, combine = TRUE)

ref_files <- list.files("/home/alan/Documents/scratchdir/downsampling/reference", full.names = TRUE)
ref_df <- read_triplets(file2triplets(ref_files), combine = TRUE)
ref_df$chain <- lapply(ref_df$file, readRDS)

edr <- mclapply(
  seq_len(nrow(ds_df)),
  function(i) {
    cat(i, "/", nrow(ds_df), "\n")
    if (isTRUE(ds_df[[i, "chains"]] == 1)) {
      return(rep(list(NULL), 3))
    }
    ind <- ref_df[["downsample_rate"]] == ds_df[[i, "downsample_rate"]]
    chain <- readRDS(ds_df[[i, "file"]])
    if (length(chain) > 1) {
      suppressMessages(
        chain <- Scalability:::combine_subposteriors(
          chain,
          subset_by = "gene",
          mc.cores = 1
        )
      )
    }
    cp <- intersect(c("nu", "s", "phi"), names(chain@parameters))
    chain@parameters[cp] <- lapply(
      chain@parameters[cp],
      function(x) {
        colnames(x) <- gsub("_Batch.*", "", colnames(x))
        x
      }
    )
    chain@parameters <- Scalability:::reorder_params(
      chain@parameters,
      gene_order = rownames(ref_df[[which(ind), "chain"]]),
      cell_order = gsub(
        "_Batch.*",
        "",
        colnames(ref_df[[which(ind), "chain"]])
      )
    )
    nsamples <- nrow(
      ref_df[[which(ind), "chain"]]@parameters[["mu"]]
    )
    chain@parameters <- lapply(
      chain@parameters,
      function(x) {
        x[seq_len(nsamples), ]
      }
    )
    suppressMessages(
      de <- BASiCS_TestDE(
        ref_df[[which(ind), "chain"]],
        chain,
        Plot = FALSE,
        PlotOffset = FALSE,
        EFDR_M = NULL,
        EFDR_D = NULL,
        EFDR_R = NULL
      )
    )
    lapply(
      de@Results,
      function(x) {
        l <- DiffRes(x)[["GeneName"]]
        if (!length(l)) NULL else l
      }
    )
  },
  mc.cores = 2
)


edr_df <- do.call(rbind, edr)
colnames(edr_df) <- c("DiffExp", "DiffDisp", "DiffResDisp")
edr_df <- as.data.frame(edr_df)
ds_df <- cbind(ds_df, edr_df)

ds_df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
  ds_df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
  function(x) sapply(x, length)
)

ds_df <- merge(ds_df, data_dims)
ds_df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
  ds_df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / ds_df[["nGenes"]],
  digits = 3
)

mdf_ds <- reshape2::melt(ds_df,
  measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_ds$variable <- gsub("pDiffExp", "mu", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffDisp", "delta", mdf_ds$variable)
mdf_ds$variable <- gsub("pDiffResDisp", "epsilon", mdf_ds$variable)
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
