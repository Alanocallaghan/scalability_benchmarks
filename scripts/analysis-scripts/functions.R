
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
  elbo <- as.data.frame(elbo, stringsAsFactors=FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}
