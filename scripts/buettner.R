library("scater")
library("BASiCS")

options(stringAsFactors = FALSE)

system("wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2805/E-MTAB-2805.processed.1.zip -O E-MTAB-2805.zip")

dir <- "downloads"
system(paste("unzip -o E-MTAB-2805.zip -d", dir))
unlink("E-MTAB-2805.zip")


erccfile <- file.path(dir, "ercc_conc.txt")
system(
  paste(
    "wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt",
    "-O", erccfile
  )
)



g2m <- read.delim(file.path(dir, "G2M_singlecells_counts.txt"))
all <- g2m
all <- all[rowSums(all[, -(1:4)]) != 0, ]
metadata <- all[, 1:4]
all <- all[, -(1:4)]
rownames(all) <- metadata[, 1]


cols <- strsplit(colnames(all), "_")
id <- sapply(cols, function(x) x[[2]])
coldata <- data.frame(id = id)
isspikes <- grepl("ERCC", metadata[, 1])

ercc_conc <- read.table(
  erccfile,
  header = TRUE,
  sep = "\t",
  fill = TRUE
)

ercc_num <- matrix(
  data = NA,
  nrow = nrow(ercc_conc),
  ncol = 1
)

ercc_num[, 1] <- (ercc_conc[, 4] * (10^(-18))) * 
  (6.0221417 * (10^23))
ercc_num <- ercc_num / 2500000
rownames(ercc_num) <- ercc_conc[, 2]

counts <- all[!isspikes, ]
libsize_drop <- isOutlier(
  Matrix::colSums(counts),
  nmads = 3,
  type = "lower",
  log = TRUE
)
feature_drop <- isOutlier(
  Matrix::colSums(counts != 0),
  nmads = 3,
  type = "lower",
  log = TRUE
)
ind_drop <- !(libsize_drop | feature_drop)
all <- all[, ind_drop]
counts <- all[!isspikes, ]

ind_expressed <- (Matrix::rowMeans(counts) >= 100 & 
  Matrix::rowMeans(counts != 0) > 0.5)
counts <- counts[ind_expressed, ]
counts <- counts[grep("^ENSMUS", rownames(counts)), ]
metadata <- metadata[ind_expressed, ]
spikes <- all[isspikes, ]
nospike <- rowSums(all[isspikes, ]) == 0
spike_input <- ercc_num[rownames(all)[isspikes], 1]
spike_input <- spike_input[!nospike]
spikes <- spikes[!nospike, ]
all <- rbind(counts, spikes)

sce <- newBASiCS_Data(
  Counts = as.matrix(all),
  SpikeInfo = data.frame(
    Name = names(spike_input), 
    Molecules = spike_input
  ),
  Tech = grepl("ERCC", rownames(all))
)
# saveRDS(sce, "datasets/buettner.rds")
colnames(sce) <- colnames(all)
buettner <- sce
# usethis::use_data(buettner, overwrite = TRUE)
saveRDS(sce, "data/buettner.rds")

