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
  lapply(files, list.files, full.names = TRUE)
}

## Convergence first
parse_elbo <- function(c) {
  c <- gsub("Chain 1:\\s+", "", c)
  normal <- grep("Drawing", c)
  abnormal <- grep("Informational", c)
  if (!length(normal)) {
    stop("Something is deeply wrong here")
  }
  ind_end <- normal - 2
  if (length(abnormal)) {
    ind_end <- abnormal - 1
  }
  elbo <- c[(grep("Begin stochastic", c) + 1):ind_end]
  ## This ain't quite right
  elbo <- gsub("MAY BE DIVERGING... INSPECT ELBO", "", elbo, fixed = TRUE)
  elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
    elbo[-c(1, grep("CONVERGED", elbo))],
    "NOTCONVERGED"
  )
  ## If both mean and median elbo converge this ends up with CONVERGED CONVERGED
  ## and therefore another column
  elbo <- gsub("(MEDIAN |MEAN )?ELBO CONVERGED", "CONVERGED", elbo)
  elbo <- gsub("(\\s+CONVERGED){2}", " CONVERGED", elbo)
  elbo <- strsplit(elbo, "\\s+")
  elbo <- do.call(rbind, elbo)
  colnames(elbo) <- elbo[1, ]
  elbo <- elbo[-1, ]
  elbo <- as.data.frame(elbo, stringsAsFactors = FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}




do_de <-function(df, ref_df, match_column, mc.cores = 2) {
  edr <- mclapply(
    seq_len(nrow(df)),
    function(i) {
      cat(i, "/", nrow(df), "\n")
      if (isTRUE(df[[i, "chains"]] == 1)) {
        return(rep(list(NULL), 3))
      }
      ind <- ref_df[[match_column]] == df[[i, match_column]]
      chain <- readRDS(df[[i, "file"]])
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
        gene_order = rownames(ref_df[[which(ind), "chain"]])
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
    mc.cores = mc.cores
  )

  edr_df <- do.call(rbind, edr)
  colnames(edr_df) <- c("DiffExp", "DiffDisp", "DiffResDisp")
  edr_df <- as.data.frame(edr_df)
  df <- cbind(df, edr_df)

  df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
    df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
    function(x) sapply(x, length)
  )

  df <- merge(df, data_dims)
  df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
    df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / df[["nGenes"]],
    digits = 3
  )
  df
}
