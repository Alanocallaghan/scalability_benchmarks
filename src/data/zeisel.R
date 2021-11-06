## Copied with modifications from MarioniLab (see https://github.com/MarioniLab/RegressionBASiCS2017/issues/1)
## https://github.com/MarioniLab/RegressionBASiCS2017/blob/fb1d833614e0469db51ec1677cabf66433f5e19e/Preprocessing/Data_preparation.R#L42

library("scater")
library("BASiCS")
library("here")

data.dir <- "downloads"

dir.create(here(data.dir), showWarnings = FALSE, recursive = TRUE)
countsfile <- file.path(data.dir, "zeisel_counts.txt")
spikesfile <- file.path(data.dir, "zeisel_spikes.txt")
erccfile <- file.path(data.dir, "ercc_conc.txt")

if (!file.exists(countsfile)) {
  download.file(
    "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt",
    destfile = countsfile
  )
}
if (!file.exists(spikesfile)) {
  download.file(
    "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt",
    destfile = spikesfile
  )
}
if (!file.exists(erccfile)) {
  download.file(
    "wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt",
    destfile = erccfile
  )
}

Zeisel <- read.delim(
  countsfile,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Meta information
Zeisel.meta <- Zeisel[1:11, ]
Zeisel.meta <- t(Zeisel.meta)
colnames(Zeisel.meta) <- Zeisel.meta[2, ]
Zeisel.meta <- Zeisel.meta[-(1:2), ]
Zeisel.meta <- as.data.frame(Zeisel.meta)

# Collect biological counts
Zeisel.bio <- Zeisel[12:nrow(Zeisel),]
rownames(Zeisel.bio) <- as.character(Zeisel.bio[, 1])
Zeisel.bio <- Zeisel.bio[, -c(1, 2)]
colnames(Zeisel.bio) <- as.character(Zeisel.meta[, 8])
Zeisel.bio <- data.matrix(Zeisel.bio[rowMeans(data.matrix(Zeisel.bio)) > 0.1, ])

# Spike-in counts
Zeisel.ERCC <- read.table(
  spikesfile,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE
)

Zeisel.ERCC <- Zeisel.ERCC[12:nrow(Zeisel.ERCC), ]
rownames(Zeisel.ERCC) <- as.character(Zeisel.ERCC[, 1])
Zeisel.ERCC <- Zeisel.ERCC[, -c(1, 2)]
colnames(Zeisel.ERCC) <- as.character(Zeisel.meta[, 8])
Zeisel.ERCC[] <- sapply(Zeisel.ERCC, as.numeric)

# Look at sizes of cell groups
table(as.character(Zeisel.meta[, 9]))


# Another cell type for downsampling - CA1
CA1.ind <- vapply(
  Zeisel.meta[, 9],
  function(x) x == "pyramidal CA1",
  logical(1)
)
CA1 <- Zeisel.bio[, CA1.ind]
Zeisel.meta <- Zeisel.meta[CA1.ind, ]
ERCC <- Zeisel.ERCC[, colnames(CA1)]

input <- rbind(CA1, ERCC)

libsize_drop <- isOutlier(
  colSums(input),
  nmads = 3,
  type = "lower",
  log = TRUE
)
feature_drop <- isOutlier(
  colSums(input != 0),
  nmads = 3,
  type = "lower",
  log = TRUE
)
input <- input[, !(libsize_drop | feature_drop)]
Zeisel.meta <- Zeisel.meta[!(libsize_drop | feature_drop), ]

bio <- input[!grepl("ERCC", rownames(input)), , drop = FALSE]
ind_expressed <- rowMeans(bio) >= 10 & rowMeans(bio != 0) > 0.5

ind_expressed_ercc <- c(
  which(ind_expressed),
  grep("ERCC", rownames(input))
)
input <- input[ind_expressed_ercc, ]

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

ERCC.num[, 1] <- (ERCC.conc[, 4] * (10^(-18))) * (6.0221417 * (10^23))
ERCC.num.final <- ERCC.num / 2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[, 2]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))], 1]
SpikeInput.1 <- data.frame(
  "Name" = names(SpikeInput),
  "Molecules" = SpikeInput,
  stringsAsFactors = FALSE
)

Data <- newBASiCS_Data(
  Counts = as.matrix(input),
  Tech = grepl("ERCC", rownames(input)),
  BatchInfo = Zeisel.meta$group,
  SpikeInfo = SpikeInput.1
)
colnames(Data) <- colnames(input)
saveRDS(Data, file = "rdata/zeisel.rds")
