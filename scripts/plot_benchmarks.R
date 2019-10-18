options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())

read_triplets <- function(files, combine = FALSE) {
  triplets <- lapply(files, list.files, full.names = TRUE)
  rows <- lapply(
    triplets, 
    function(x) {
      row <- readRDS(x[[2]])
      row <- as.data.frame(row)
      row[["time"]] <- readRDS(x[[3]])[["elapsed"]]
      row
    }
  )
  df <- do.call(rbind, rows)
  df[["file"]] <- lapply(triplets, function(x) x[[1]])
  # if (combine) {
  #   df[["chain"]] <- lapply(
  #     df[["file"]], 
  #     function(x) {
  #       c <- readRDS(x)
  #       if (length(c) > 1) {
  #         combine_subposteriors(c, subset_by = "gene")
  #       } else {
  #         c
  #       }
  #     }
  #   )
  # } else {
  #   df[["chain"]] <- lapply(df[["file"]], readRDS)
  # }
  df
}


advi_files <- list.files("/home/alan/Documents/scratchdir/advi", full.names = TRUE)
advi_df <- read_triplets(advi_files)



dc_files <- list.files("/home/alan/Documents/scratchdir/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(dc_files, combine = TRUE)

df <- rbind(advi_df, dc_df)

datasets <- unique(df[["data"]])
data_dims <- vapply(
  datasets,
  function(x) {
    suppressMessages(
      dim(
        readRDS(paste0("data/", x, ".rds"))
      )
    )
  },
  FUN.VALUE = numeric(2)
)
data_dims <- as.data.frame(t(data_dims))
colnames(data_dims) <- c("nGenes", "nCells")
data_dims[["data"]] <- datasets
df <- merge(df, data_dims)

source("scripts/time_plot.R")
references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

source("scripts/de_on_table.R")
