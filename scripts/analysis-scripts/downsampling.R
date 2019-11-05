options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())
source(here("scripts/analysis-scripts/functions.R"))

ds_files <- list.files("/home/alan/Documents/scratchdir/downsampling", full.names = TRUE)
ds_df <- read_triplets(file2triplets(ds_files), combine = TRUE)

edr <- mclapply(
  seq_len(nrow(ds_df)),
  function(i) {
    cat(i, "/", nrow(ds_df), "\n")
    if (isTRUE(ds_df[[i, "chains"]] == 1)) {
      return(rep(list(NULL), 3))
    }
    ind <- references[["data"]] == ds_df[[i, "data"]]
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
      gene_order = rownames(references[[which(ind), "chain"]]),
      cell_order = gsub(
        "_Batch.*",
        "",
        colnames(references[[which(ind), "chain"]])
      )
    )
    nsamples <- nrow(
      references[[which(ind), "chain"]]@parameters[["mu"]]
    )
    chain@parameters <- lapply(
      chain@parameters,
      function(x) {
        x[seq_len(nsamples), ]
      }
    )
    suppressMessages(
      de <- BASiCS_TestDE(
        references[[which(ind), "chain"]],
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

ggplot(mdf_ds, aes(x = downsample_rate, y = value, color = variable)) +
  geom_point() +
  scale_color_brewer(name = "Parameter", palette = "Set1") +
  scale_x_continuous(label = scales::percent) +
  labs(x = "Downsample rate", y = "Portion of genes perturbed")
