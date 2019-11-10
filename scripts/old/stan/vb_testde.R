library("here")
library("devtools")

load_all("../BASiCS")
load_all(here())
options(mc.cores = 4)


d <- zd()
# d <- makeExampleBASiCS_Data()

spikes <- d[isSpike(d), ]
counts <- d[!isSpike(d), ]
spikes <- assay(spikes)
counts <- assay(counts)

L <- 12


start <- HiddenBASiCS_MCMC_Start(d, 
  eta = 5, 
  m = rep(0, L), 
  V = diag(L), 
  a.sigma2 = 2, 
  b.sigma2 = 2, 
  WithSpikes = TRUE)

ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)
sdata <- list(
  G = nrow(counts), 
  C = ncol(counts),
  SG = nrow(spikes),
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
  L = 12,
  vbeta = diag(L),
  rbf_variance = 1.2,
  eta = 5,
  astwo = 2,
  bstwo = 2
)


set.seed(42)
suppressMessages(
  tmp <- capture.output(
    b <- BASiCS_MCMC(
      d, 
      N = 20000, 
      Thin = 10, 
      Burn = 10000,
      Regression = TRUE,
      WithSpikes = TRUE,
      FixML = TRUE,
      ml = ml,
      PrintProgress = FALSE)
  )
)

 
vb_objects <- mclapply(benchmarkSeeds(), function(seed) {
  set.seed(seed)
  c <- capture.output(v <- vb(stanmodels$basics_regression, iter = 100000, data = sdata))
  list(c, v)
}, mc.cores = 6)

# b_nr <- BASiCS_MCMC(
#   d, 
#   N = 20000, 
#   Thin = 10, 
#   Burn = 10000,
#   Regression = TRUE,
#   WithSpikes = TRUE,
#   FixML = TRUE,
#   ml = ml,
#   PrintProgress = FALSE)

# vb_objects_nr <- lapply(benchmarkSeeds(), function(seed) {
#   set.seed(seed)
#   v <- vb(stanmodels$basics, iter = 100000, data = sdata)
#   v
# })

saveRDS(b, "outputs/basics.rds")
saveRDS(vb_objects, "outputs/vb_objects.rds")


# saveRDS(b, "outputs/basics_nr.rds")
# saveRDS(vb_objects, "outputs/vb_objects_nr.rds")

# saveRDS(resdf, "outputs/vb_testde.rds")





# vb_b_objects <- lapply(vb_objects, function(x) {
#   xe <- extract(x)
#   parameters <- list(
#     mu = xe$mu,
#     delta = xe$delta,
#     epsilon = xe$epsilon,
#     s = xe$s,
#     nu = xe$nu,
#     theta = as.matrix(xe$theta),
#     phi = xe$phi
#   )
#   gp <- c("mu", "delta", "epsilon")
#   cp <- c("s", "nu", "phi")
#   parameters[gp] <- lapply(gp, function(x) {
#     colnames(parameters[[x]]) <- colnames(b@parameters[["mu"]])
#     parameters[[x]]
#   })
#   parameters[cp] <- lapply(cp, function(x) {
#     colnames(parameters[[x]]) <- colnames(b@parameters[["nu"]])
#     parameters[[x]]
#   })

#   new("BASiCS_Chain", parameters = parameters)
# })


# res <- lapply(seq_along(vb_b_objects), function(i) {
#   vb <- vb_b_objects[[i]]
#   de <- BASiCS_TestDE(
#     vb, 
#     b, 
#     GroupLabel1 = "VB", 
#     GroupLabel2 = "MCMC", 
#     Plot = FALSE,
#     PlotOffset = FALSE
#   )
#   diffexp <- function(foo) {
#     !foo %in% c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")
#   }  
#   data.frame(
#     Ind = i,
#     "DiffMean" = sum(diffexp(de$TableMean$ResultDiffMean)),
#     "DiffDisp" = sum(diffexp(de$TableMean$ResultDiffMean)),
#     "DiffResDisp" = sum(diffexp(de$TableResDisp$ResultDiffResDisp))
#   )
# })

# resdf <- do.call(rbind, res)

# colnames(resdf) <- c(
#   "Instance", 
#   "Diff. expressed", 
#   "Diff. variability", 
#   "Diff. residual variability"
# )
# resdf <- rbind(resdf, c("Mean", round(colMeans(resdf[, 2:4]), digits = 2)))



# # max fold change vs chains
# fc <- function(chain, param) {
#   colMedians(chain@parameters[[param]]) / 
#     colMedians(c1@parameters[[param]])
# }



# sapply(vb_b_objects, function(x) max(fc(x, "mu")))
# sapply(vb_b_objects, function(x) max(fc(x, "delta")))
# sapply(vb_b_objects, function(x) max(fc(x, "epsilon")))



# print(xtable::xtable(resdf, 
#   label = "tab:DiffResVB",
#   caption = "Table showing detected differentially expressed and differentially variable genes in different runs of ADVI.",
#   ), 
#   include.rownames = FALSE, 
#   file = "/home/alan/Documents/github/latex-documents/tables/scalability/vb_TestDE.tex")

