library("BASiCS")
library("reshape2")
library("xtable")

out <- readRDS("outputs/hmc_vs_amwg.rds")

hmc_ess <- sapply(c("mu", "delta", "epsilon"),
    function(p) mean(BASiCS_EffectiveSize(out$chains$hmc, p) / out$time$hmc[["elapsed"]])
)
amwg_ess <- sapply(c("mu", "delta", "epsilon"),
    function(p) mean(BASiCS_EffectiveSize(out$chains$amwg, p) / out$time$amwg[["elapsed"]])
)

data <- cbind(
    HMC = hmc_ess,
    aMwG = amwg_ess
)
data <- melt(data)
colnames(data) <- c("Parameter", "Method", "ESS/s")
data <- data[, c("Method", "Parameter", "ESS/s")]
table <- data[c(1, 4, 2, 5, 3, 6), ]
fdf <- format(table, digits = 1)

tab <- xtable(
    fdf,
    align = "rr|rr",
    caption = "Average effective sample size per second for HMC and AMWG samplers using \\cite{Tung2017} data.",
    label = "tab:HMCvsAMWG"
)
print(tab,
    include.rownames = FALSE,
    hline.after = 0,
    sanitize.text.function = function(x) {x},
    file = "tables/hmc-comparison.tex"
)
