suppressPackageStartupMessages({
  library("rstan")
  library("here")
  library("devtools")
  load_all(here("../BASiCS"))
  load_all(here())
})
source(here("data-raw/stan/functions.R"))


tungData <- function(cell_ind = 1) {

  reads <- read.delim(
    file.path(system.file("datasets/tung_molecules-filter.txt", package = "Scalability")),
    stringsAsFactors = FALSE
  )
  rownames(reads) <- reads[[1]]
  reads <- reads[, -1]
  reads <- as.matrix(reads)

  cell_lines <- gsub("(.*)\\.r\\d\\..*", "\\1", colnames(reads))
  ind_cell_one <- cell_lines == sort(unique(cell_lines))[[cell_ind]]

  reads <- reads[, ind_cell_one]

  reads <- reads[rowSums(reads) > 0, ]
  ind_expressed <- rowMeans(reads) >= 1
  reads <- reads[ind_expressed, ]
  reads <- reads[, colMeans(reads !=0) > 0.65]

  libsize_drop <- isOutlier(colSums(reads), nmads = 3, type = "lower", log = TRUE)
  feature_drop <- isOutlier(colSums(reads != 0), nmads = 3, type = "lower", log = TRUE)

  reads <- reads[, !(libsize_drop | feature_drop)]

  spikes <- rownames(reads)[grep("ERC", rownames(reads))]


  ERCC.conc <- read.table(
    paste0(
      file.path(system.file("datasets", package = "Scalability"), "ercc_conc.txt")
    ),
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

  SpikeInput <- ERCC.num.final[rownames(reads)[grepl("ERCC", rownames(reads))], 1]
  SpikeInput.1 <- data.frame(
    "Name" = names(SpikeInput),
    "Molecules" = SpikeInput,
    stringsAsFactors = FALSE)

  batches <- gsub(".*\\.(r\\d)\\..*", "\\1", colnames(reads))

  bd <- newBASiCS_Data(
    Counts = reads,
    BatchInfo = batches,
    Tech = rownames(reads) %in% spikes,
    SpikeInfo = SpikeInput.1
  )
  bd
}

d <- tungData()

spikes <- d[isSpike(d), ]
counts <- d[!isSpike(d), ]
spikes <- assay(spikes)
counts <- assay(counts)

l <- 12


start <- BASiCS:::HiddenBASiCS_MCMC_Start(
  d, 
  eta = 5, 
  m = rep(0, l), 
  V = diag(l), 
  a.sigma2 = 2, 
  b.sigma2 = 2, 
  WithSpikes = TRUE)

ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)
sdata <- list(
  g = nrow(counts), 
  c = ncol(counts),
  sg = nrow(spikes),
  counts = counts, 
  spikes = spikes,
  spike_levels = metadata(d)$SpikeInput,
  as = 1,
  bs = 1,
  atheta = 1,
  btheta = 1,
  smu = sqrt(0.5),
  sdelta = sqrt(0.5),
  aphi = rep(1, ncol(counts)),
  mbeta = rep(0, L),
  ml = ml,
  l = 12,
  vbeta = diag(L),
  rbf_variance = 1.2,
  eta = 5,
  astwo = 2,
  bstwo = 2
)


# tm <- system.time({
#   set.seed(42)
#   m <- sampling(
#     stanmodels$basics, 
#     data = sdata)  
# })
# saveRDS(m, file = here("outputs/hmc_tung.rds"))

tmr <- system.time({
  set.seed(42)
  m_reg <- sampling(
    stanmodels$basics_regression, 
    data = sdata)  
})

saveRDS(m, file = here("outputs/hmc_tung_reg.rds"))
rm(m)
rm(m_reg)


tv <- system.time({
  set.seed(42)
  v <- vb(stanmodels$basics, iter = 1000000, data = sdata)    
})
saveRDS(v, file = here("outputs/vb_tung.rds"))


tvr <- system.time({
  set.seed(42)
  v_reg <- vb(stanmodels$basics_regression, iter = 1000000, data = sdata)    
})
saveRDS(v_reg, file = here("outputs/vb_tung_reg.rds"))
rm(v)
rm(v_reg)


tb <- system.time({
  set.seed(42)
  b <- BASiCS_MCMC(
    d,
    ml  = ml,
    N = 20000, 
    Thin = 10, 
    Burn = 10000, 
    Regression = FALSE, 
    WithSpikes = TRUE)
})

saveRDS(b, file = here("outputs/basics_tung.rds"))

tbr <- system.time({
  set.seed(42)
  b_reg <- BASiCS_MCMC(
    d,
    ml  = ml,
    N = 20000, 
    Thin = 10, 
    Burn = 10000, 
    Regression = TRUE, 
    WithSpikes = TRUE)
})

saveRDS(b_reg, file = here("outputs/basics_tung_reg.rds"))

rm(b)
rm(b_reg)
