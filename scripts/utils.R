read_triplets <- function(files) {
  triplets <- lapply(
    files,
    function(file) {
      con <- file(file)
      text <- readLines(con, warn = FALSE)
      close(con)
      text <- gsub("(\\[|\\])", "", text)
      text <- strsplit(text, ", ")[[1]]
    }
  )
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
  # df[["chain"]] <- lapply(df[["file"]], readRDS)
  df
}
