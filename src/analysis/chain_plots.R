dir.create("figs/ess", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/hpd", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/fit", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/de", showWarnings = FALSE, recursive = TRUE)



## Mean-variance curves and DE plots for worst for each data
maxdfm <- df %>%
  filter(!is.na(chains)) %>%
  group_by(data) %>%
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)

maxdfe <- df %>%
  filter(!is.na(chains)) %>%
  group_by(data) %>%
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)

cat("ESS\n")
ess_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_ess_plot,
  maxdf = maxdfe,
  references = references
)

cat("Fits\n")
fit_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_fit_plot,
  maxdf = maxdfe,
  references = references
)

cat("DE\n")
## MA plot of worst case for worst for each data
de_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_de_plot,
  maxdf = maxdfe,
  references = references
)

cat("HPD\n")
hpd_plots <- lapply(
  seq_len(nrow(maxdfe)),
  do_hpd_plots,
  maxdf = maxdfe,
  references = references
)

## Posterior variance plots
maxdfadvim <- df %>% filter(is.na(chains)) %>%
  group_by(data) %>%
  top_n(n = 1, wt = nDiffExp) %>%
  distinct(data, .keep_all = TRUE)

maxdfadvie <- df %>% filter(is.na(chains)) %>%
  group_by(data) %>%
  top_n(n = 1, wt = nDiffResDisp) %>%
  distinct(data, .keep_all = TRUE)


cat("Fit ADVI\n")
fit_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_fit_plot,
  maxdf = maxdfadvie,
  references = references
)

cat("DE ADVI\n")
## MA plot of worst case for worst for each data
de_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_de_plot,
  maxdf = maxdfadvie,
  references = references
)

cat("HPD ADVI\n")
hpd_plots_advi <- lapply(
  seq_len(nrow(maxdfadvie)),
  do_hpd_plots,
  maxdf = maxdfadvie,
  references = references
)
