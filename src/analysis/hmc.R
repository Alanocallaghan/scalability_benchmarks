library("BASiCS")
library("reshape2")
library("xtable")

out <- readRDS("outputs/hmc_vs_amwg.rds")


hmc_ess <- lapply(c("mu", "delta", "epsilon"),
    function(p) mean(BASiCS_EffectiveSize(out$chains$hmc, p) / out$time$hmc[["elapsed"]])
)
amwg_ess <- lapply(c("mu", "delta", "epsilon"),
    function(p) mean(BASiCS_EffectiveSize(out$chains$amwg, p) / out$time$amwg[["elapsed"]])
)

data <- rbind(
    do.call(cbind, hmc_ess),
    do.call(cbind, amwg_ess)
)
colnames(data) <- c("mu", "delta", "epsilon")
table <- data.frame(
    Method = c("AMWG", "HMC"),
    data
)
mdf <- melt(table, id.var = "Method")
colnames(mdf) <- c("Method", "Parameter", "ESS/s")
fdf <- format(mdf, digits = 1)

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
