library("scater")
library("BASiCS")
library("here")



data.dir <- "downloads"
tungfile <- file.path(data.dir, "tung_molecules-filter.txt")
erccfile <- file.path(data.dir, "ercc_conc.txt")
system(
  paste(
    "wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt",
    "-O", erccfile
  )
)

system(
  paste(
    "wget", 
    "https://github.com/Alanocallaghan/singleCellSeq/raw/master/data/molecules-filter.txt",
    "-O", tungfile
  )
)

reads <- read.delim(
  tungfile,
  stringsAsFactors = FALSE
)
rownames(reads) <- reads[[1]]
reads <- reads[, -1]
reads <- as.matrix(reads)

cell_lines <- gsub("(.*)\\.r\\d\\..*", "\\1", colnames(reads))
cell_ind <- 1
ind_cell_one <- cell_lines == sort(unique(cell_lines))[[cell_ind]]

reads <- reads[, ind_cell_one]

reads <- reads[rowSums(reads) > 0, ]
ind_expressed <- rowMeans(reads) >= 1 & rowMeans(reads != 0) > 0.5
reads <- reads[ind_expressed, ]
reads <- reads[, colMeans(reads != 0) > 0.65]

libsize_drop <- isOutlier(colSums(reads), nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(colSums(reads != 0), nmads = 3, type = "lower", log = TRUE)

reads <- reads[, !(libsize_drop | feature_drop)]

spikes <- rownames(reads)[grep("ERC", rownames(reads))]


ERCC.conc <- read.table(
  erccfile,
  header = TRUE,
  sep = "\t",
  fill = TRUE
)

ERCC.num <- matrix(
  data = NA,
  nrow = nrow(ERCC.conc),
  ncol = 1
)

ERCC.num[, 1] <- (ERCC.conc[, 4] * (10^(-18))) * 
  (6.0221417 * (10^23))
ERCC.num.final <- ERCC.num / 2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[, 2]

SpikeInput <- ERCC.num.final[rownames(reads)[grepl("ERCC", rownames(reads))], 1]
SpikeInput.1 <- data.frame(
  "Name" = names(SpikeInput),
  "Molecules" = SpikeInput,
  stringsAsFactors = FALSE
)

batches <- gsub(".*\\.(r\\d)\\..*", "\\1", colnames(reads))

bd <- newBASiCS_Data(
  Counts = reads,
  BatchInfo = batches,
  Tech = rownames(reads) %in% spikes,
  SpikeInfo = SpikeInput.1
)
colnames(bd) <- colnames(reads)
saveRDS(bd, file = "data/tung.rds")
# usethis::use_data(tung, overwrite = TRUE)
# saveRDS(bd, "datasets/tung.rds")
