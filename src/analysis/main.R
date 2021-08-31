library("here")

source(here("src/analysis/preamble.R"))

advi_files <- list.files("outputs/advi", full.names = TRUE)
advi_triplets <- file2triplets(advi_files)
advi_elbo <- lapply(advi_triplets, function(x) readRDS(x[[3]]))
advi_triplets <- lapply(advi_triplets, function(x) x[-3])
advi_df <- read_triplets(advi_triplets)




dc_files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
dc_df <- read_triplets(file2triplets(dc_files), combine = TRUE)


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
