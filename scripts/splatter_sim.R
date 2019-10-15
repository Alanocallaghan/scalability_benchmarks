library("splatter")
devtools::load_all()

d <- zeiselData()
params <- splatEstimate(d)
sim <- splatSimulate(params, batchCells = c(5000, 5000))
isSpike(d, "ERCC")

ercc_conc <- read.table(
  paste0(
    file.path(system.file("datasets", package = "Scalability"), "ercc_conc.txt")
  ),
  header = TRUE,
  sep = "\t",
  fill = TRUE
)

ercc_num <- matrix(
  data = NA,
  nrow = nrow(ercc_conc),
  ncol = 1
)

ercc_num[, 1] <- (ercc_conc[, 4] * (10^(-18))) * (6.0221417 * (10^23))
ercc_num <- ercc_num / 2500000
rownames(ercc_num) <- ercc_conc[, 2]

spike_input <- ercc_num[rownames(d)[grepl("ERCC", rownames(d))], 1]
metadata(sim)$SpikeInput <- spike_input
isSpike(sim, "ERCC") <- grepl("ERCC", rownames(d))
rownames(sim)[isSpike(d)] <- rownames(d)[isSpike(d)]
colData(sim)$BatchInfo <- colData(sim)$Batch
splatter <- sim
usethis::use_data(splatter)
# saveRDS(sim, file = "datasets/splatter_sim.rds")
