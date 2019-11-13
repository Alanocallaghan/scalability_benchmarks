options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())
source(here("scripts/analysis-scripts/functions.R"))

rm_files <- list.files("outputs/removing_cells/removing/", full.names = TRUE)
rmt <- file2triplets(rm_files)
rmt <- rmt[as.logical(sapply(rmt, length))]
rm_df <- read_triplets(rmt, combine = TRUE)

rm_ref_files <- list.files("outputs/removing_cells/reference", full.names = TRUE)
rm_ref_df <- read_triplets(file2triplets(rm_ref_files), combine = TRUE)
rm_ref_df$chain <- lapply(rm_ref_df$file, readRDS)

edr <- mclapply(
  seq_len(nrow(rm_df)),
  function(i) {
    cat(i, "/", nrow(rm_df), "\n")
    if (isTRUE(rm_df[[i, "chains"]] == 1)) {
      return(rep(list(NULL), 3))
    }
    ind <- rm_ref_df[["proportion_retained"]] == rm_df[[i, "proportion_retained"]]
    chain <- readRDS(rm_df[[i, "file"]])
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
      gene_order = rownames(rm_ref_df[[which(ind), "chain"]])
    )
    nsamples <- nrow(
      rm_ref_df[[which(ind), "chain"]]@parameters[["mu"]]
    )
    chain@parameters <- lapply(
      chain@parameters,
      function(x) {
        x[seq_len(nsamples), ]
      }
    )
    suppressMessages(
      de <- BASiCS_TestDE(
        rm_ref_df[[which(ind), "chain"]],
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
rm_df <- cbind(rm_df, edr_df)

rm_df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
  rm_df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
  function(x) sapply(x, length)
)

rm_df <- merge(rm_df, data_dims)
rm_df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
  rm_df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / rm_df[["nGenes"]],
  digits = 3
)

mdf_rm <- reshape2::melt(rm_df,
  measure.vars = c("pDiffExp", "pDiffDisp", "pDiffResDisp")
)
mdf_rm$variable <- gsub("pDiffExp", "mu", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffDisp", "delta", mdf_rm$variable)
mdf_rm$variable <- gsub("pDiffResDisp", "epsilon", mdf_rm$variable)
mdf_rm$cells_retained <- mdf_rm$proportion_retained * mdf_rm$nCells

mdf_rm$proportion_retained <- factor(
  paste(mdf_rm$proportion_retained * 100, "%"),
  levels = paste(sort(unique(rm_df$proportion_retained), decreasing = TRUE) * 100, "%")
)


ggplot(mdf_rm, aes(x = cells_retained, y = value, color = variable)) +
  geom_quasirandom(dodge.width = 30, size = 0.3) +
  scale_color_brewer(name = "Parameter", palette = "Set1") +
  scale_x_reverse() +
  scale_y_continuous(label = scales::percent) +
  labs(x = "nCells", y = "Portion of genes perturbed")

ggsave("figs/removing_cells.pdf", width = 6, height = 4)
