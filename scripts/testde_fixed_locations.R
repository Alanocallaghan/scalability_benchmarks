library("here")
library("devtools")

load_all("../BASiCS")
load_all(here())



test_fixed <- function(d) {
  L <- 12
  start <- HiddenBASiCS_MCMC_Start(d, 
    eta = 5, 
    m = rep(0, L), 
    V = diag(L), 
    a.sigma2 = 2, 
    b.sigma2 = 2, 
    WithSpikes = TRUE)

  ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)

  b1 <- BASiCS_MCMC(
    d,
    N = 20000, Thin = 10, Burn = 10000,
    Regression = TRUE, 
    WithSpikes = TRUE,
    PrintProgress = FALSE,
    ml = ml)


  b2 <- BASiCS_MCMC(
    d,
    N = 20000, Thin = 10, Burn = 10000,
    Regression = TRUE, 
    PrintProgress = FALSE,
    WithSpikes = TRUE)


  BASiCS_TestDE(b1, b2, 
    GroupLabel1 = "Fixed", GroupLabel2 = "Not fixed",
    Plot = FALSE, PlotOffset = FALSE)
}


b <- readRDS("outputs/basics_zeisel.rds")

d <- zd()
L <- 12
start <- HiddenBASiCS_MCMC_Start(d, 
  eta = 5, 
  m = rep(0, L), 
  V = diag(L), 
  a.sigma2 = 2, 
  b.sigma2 = 2, 
  WithSpikes = TRUE)

ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)

bf <- readRDS("outputs/basics_zeisel_fixed.rds")
# bf <- BASiCS_MCMC(
#   d,
#   N = 20000, Thin = 10, Burn = 10000,
#   Regression = TRUE, 
#   WithSpikes = TRUE,
#   PrintProgress = FALSE,
#   ml = ml)





d <- zd()
L <- 12
start <- HiddenBASiCS_MCMC_Start(d, 
  eta = 5, 
  m = rep(0, L), 
  V = diag(L), 
  a.sigma2 = 2, 
  b.sigma2 = 2, 
  WithSpikes = TRUE)

ml <- BASiCS:::HiddenFindRBFLocations(start$mu0, L)

bf <- BASiCS_MCMC(
  d,
  N = 20000, Thin = 10, Burn = 10000,
  Regression = TRUE, 
  WithSpikes = TRUE,
  PrintProgress = FALSE,
  ml = ml)

saveRDS(b1, "outputs/basics_zeisel_fixed.rds")


tung_res <- test_fixed(td())
zeis_res <- test_fixed(zd())

saveRDS(tung_res, "outputs/tung_fixed_de.rds")
saveRDS(zeis_res, "outputs/zeisel_fixed_de.rds")


diffexp <- function(foo) {
  !foo %in% c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")
}  

tung_tab <- data.frame(
  dataset = "Tung",
  "DiffMean" = sum(diffexp(tung_res$TableMean$ResultDiffMean)),
  "DiffDisp" = sum(diffexp(tung_res$TableDisp$ResultDiffDisp)),
  "DiffResDisp" = sum(diffexp(tung_res$TableResDisp$ResultDiffResDisp))
)

zeis_tab <- data.frame(
  dataset = "Zeisel",
  "DiffMean" = sum(diffexp(zeis_res$TableMean$ResultDiffMean)),
  "DiffDisp" = sum(diffexp(zeis_res$TableDisp$ResultDiffDisp)),
  "DiffResDisp" = sum(diffexp(zeis_res$TableResDisp$ResultDiffResDisp))
)

# tab <- rbind(zeis_tab, tung_tab)
tab <- zeis_tab
colnames(tab) <- c("Dataset", "Diff. mean", "Diff. variability", "Diff. residual variability")
tab <- tab[, -1]

print(xtable::xtable(tab, 
  label = "table::FixedLocationDiff",
  caption = "Number of genes highlighted as having differential expression and 
    differential (residual) over-dispersion using BASiCS with adaptive GRBF 
    locations, relative to BASiCS using fixed GRBF locations.",
  ), 
  include.rownames = FALSE, 
  file = "/home/alan/Documents/github/latex-documents/tables/scalability/fixedLocationDiff.tex"
)







