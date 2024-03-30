library("here")

source(here("src/analysis/preamble.R"))
source(here("src/analysis/read_chains.R"))

sourceme(here("src/analysis/de_on_table.R"))
sourceme(here("src/analysis/chain_plots.R"))

## these can be run separately...
# sourceme(here("src/analysis/diagnostics.R"))
# sourceme(here("src/analysis/normalisation_comparison.R"))
# sourceme(here("src/analysis/point_estimates.R"))
# sourceme(here("src/analysis/hpd.R"))

## not doing this any more...
## sourceme(here("src/analysis/identifiability.R"))

save.image("main_done.RData")
