library("devtools")
load_all("../BASiCS")
load_all()


d <- makeExampleBASiCS_Data(WithSpikes = FALSE)

d <- williamsData()

m <- BASiCS_MCMC(d,
    N = 20000,
    Thin = 10,
    Burn = 10000,
    Regression = TRUE,
    WithSpikes = FALSE,
    PrintProgress = TRUE
)

system.time(sv <- BASiCS_stan(d, regression = TRUE, with_spikes = FALSE))
sm <- BASiCS_stan(d, method = "sampling", regression = TRUE, with_spikes = FALSE, verbose = TRUE)

saveRDS(sv, "advi_williams.rds")
saveRDS(sm, "hmc_williams.rds")
saveRDS(m, "basics_williams.rds")



svb <- stan2basics(sv, 
    gene_names = colnames(m@parameters[["mu"]]),
    cell_names = colnames(m@parameters[["nu"]])
)
smb <- stan2basics(sv, 
    gene_names = colnames(m@parameters[["mu"]]),
    cell_names = colnames(m@parameters[["nu"]])
)

svb <- offset_correct(m, svb)
smb <- offset_correct(m, smb)


comp_plot(m, svb, smb, names = c("AMWG", "ADVI", "HMC"))



comp_plot(m, svb, names = c("AMWG", "ADVI"))
