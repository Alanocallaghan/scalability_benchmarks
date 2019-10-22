options(stringsAsFactors = FALSE)
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("BASiCS")
library("Scalability")

theme_set(theme_bw())

read_triplets <- function(triplets, combine = FALSE) {
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
  df
}

file2triplets <- function(files) {
  triplets <- lapply(files, list.files, full.names = TRUE)
}

## Convergence first
parse_elbo <- function(c) {
  c <- gsub("Chain 1:\\s+", "", c)
  elbo <- c[(grep("Begin stochastic", c) + 1):(grep("Drawing", c) - 2)]
  elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
    elbo[-c(1, grep("CONVERGED", elbo))],
    "NOTCONVERGED")
  elbo <- gsub("(MEDIAN )?ELBO CONVERGED", "CONVERGED", elbo)
  elbo <- strsplit(elbo, "\\s+")
  elbo <- do.call(rbind, elbo)
  colnames(elbo) <- elbo[1, ]
  elbo <- elbo[-1, ]
  elbo <- as.data.frame(elbo, stringsAsFactors=FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}


advi_files <- list.files("/home/alan/Documents/scratchdir/advi", full.names = TRUE)
advi_triplets <- file2triplets(advi_files)
advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
advi_triplets <- lapply(advi_triplets, function(x) x[-3])
advi_df <- read_triplets(advi_triplets)


dc_files <- list.files("/home/alan/Documents/scratchdir/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)

file_df <- rbind(advi_df, dc_df)

datasets <- unique(file_df[["data"]])
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
df <- merge(file_df, data_dims)

source(here("scripts/analysis-scripts/time_plot.R"))
references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

source(here("scripts/analysis-scripts/de_on_table.R"))
