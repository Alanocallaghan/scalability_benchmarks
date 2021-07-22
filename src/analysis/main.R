library("here")

source(here("src/analysis/preamble.R"))

file_df <- rbind(advi_df, dc_df)
# file_df <- dc_df
df <- merge(file_df, data_dims)

references <- df[which(df[["chains"]] == 1), ]
references[["chain"]] <- lapply(references[["file"]], readRDS)

sourceme(here("src/analysis/de_on_table.R"))
sourceme(here("src/analysis/chain_plots.R"))

sourceme(here("src/analysis/normalisation_comparison.R"))

## not doing this any more...
## sourceme(here("src/analysis/identifiability.R"))
