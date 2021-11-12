###############################################################################
##
## Run TestDE
##
###############################################################################

df <- do_de(df, ref_df = references, match_column = "data", data_dims,
  mc.cores = 2)


###############################################################################
##
## DE numbers
##
###############################################################################

sdf  <- df[,
  c("data", "chains", "pDiffExp", "pDiffDisp", "pDiffResDisp")
]
mdf <- reshape2::melt(sdf, id.vars = c("data", "chains"))
mdf$chains <- as.numeric(mdf$chains)
mdf$variable <- gsub("^pDiffExp$", "mu", mdf$variable)
mdf$variable <- gsub("^pDiffDisp$", "delta", mdf$variable)
mdf$variable <- gsub("^pDiffResDisp$", "epsilon", mdf$variable)
mdf$variable <- factor(
  mdf$variable,
  levels = c("mu", "delta", "epsilon")
)

mdf$data <- sub(
  "([[:alpha:]])", "\\U\\1",
  mdf$data,
  perl = TRUE
)
mdf$data <- sub(
  "Pbmc",
  "10x PBMC",
  mdf$data
)

# mdf <- mdf %>% filter(variable %in% c("mu", "epsilon"))

advi_mdf <- mdf %>%
  filter(is.na(chains)) %>%
  group_by(data, variable, chains) %>%
  summarize(value = median(value))


g <- ggplot(mdf[!(is.na(mdf$chains) | mdf$chains == 1), ],
  aes(
    x = factor(chains),
    y = value,
    # group = data,
    color = variable
  )
) +
  geom_quasirandom(
    groupOnX = TRUE,
    dodge.width = 0.5,
    # position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
    size = 0.5
  ) +
  # geom_point(
  #   position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
  #   size = 0.5
  # ) +
  # geom_violin() +
  geom_hline(
    data = advi_mdf,
    alpha = 0.5,
    aes(
      yintercept = value,
      color = variable,
      linetype = "ADVI"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    labels = "Mean\nADVI results",
    values = 2
  ) +
  facet_wrap(~data, nrow = 2, ncol = 2) +
  scale_x_discrete(name = "Partitions") +
  scale_y_continuous(
    name = "Portion of genes perturbed",
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  theme(text = element_text(size = 18)) +
  scale_color_brewer(name = "Parameter", palette = "Set1")

ggsave(here("figs/diffexp_plot.pdf"), width = 12, height = 8)

###############################################################################
##
## DE overlaps
##
###############################################################################

md <- df %>% group_by(chains, data, by)


count_instances <- function(x) {
  union <- Reduce(union, x)
  counts <- sapply(
    union,
    function(y) {
      sum(sapply(x, function(z) y %in% z))
    }
  )
  list(setNames(counts, union))
}


mds <- md %>% dplyr::summarise(
  .groups = "drop_last",
  mu = count_instances(DiffExp),
  delta = count_instances(DiffDisp),
  epsilon = count_instances(DiffResDisp)
)




overlap_df <- function(mds, i, var) {
  o <- mds[[i, var]]
  if (length(o) > 0 && !all(is.na(names(o))) && length(o)[[1]] > 0) {
    data.frame(
      by = mds[[i, "by"]],
      chains = mds[[i, "chains"]],
      data = mds[[i, "data"]],
      var = var,
      n = o
    )
  } else {
    matrix(
      nrow = 0,
      ncol = 5,
      dimnames = list(NULL, c("by", "chains", "data", "var", "n"))
    )
  }
}

all_overlap_df <- lapply(seq_len(nrow(mds)), function(i) {
  d1 <- overlap_df(mds, i, "mu")
  d2 <- overlap_df(mds, i, "delta")
  d3 <- overlap_df(mds, i, "epsilon")
  as.data.frame(rbind(d1, d2, d3))
})

all_overlap_df <- do.call(rbind, all_overlap_df)
all_overlap_df$Gene <- rownames(all_overlap_df)
all_overlap_df$n <- factor(all_overlap_df$n)
all_overlap_df <- merge(data_dims, all_overlap_df, by = "data")
all_overlap_df$data <- sub(
  "([[:alpha:]])([[:alpha:]]+)",
  "\\U\\1\\L\\2",
  all_overlap_df$data,
  perl = TRUE
)
all_overlap_df$data <- gsub("Pbmc", "10X PBMC", all_overlap_df$data)
if (nrow(all_overlap_df)) {
  all_overlap_df$chains <- paste(all_overlap_df$chains, "chains")
  levs <- c(
    paste(
      sort(as.numeric(unique(all_overlap_df$chains))), "chains"),
      "ADVI"
  )
} else {
  levs <- c("2 chains", "ADVI")
}
all_overlap_df$chains <- factor(
  all_overlap_df$chains,
  levels = levs
)
all_overlap_df$data <- factor(all_overlap_df$data,
  levels = unique(mdf$data)
)
all_overlap_df$chains[is.na(all_overlap_df$chains)] <- "ADVI"


# all_overlap_df <- all_overlap_df[!is.na(all_overlap_df$chains), ]


count_df <- all_overlap_df %>%
  group_by(n, data, chains, var, nGenes) %>%
  summarise(count = n())

## Proportions because of differing number of genes
g <- ggplot(
  count_df,
    aes(
      x = n,
      fill = var,
      y = count / nGenes
    )
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(preserve = "single")
  ) +
  labs(
    x = "Number of times gene was identified",
    y = "Portion of genes"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  scale_fill_brewer(
    palette = "Set1", name = "Parameter",
    limits = c("mu", "delta", "epsilon"),
    labels = c(bquote(mu[i]), bquote(delta[i]), bquote(epsilon[i]))) +
  facet_grid(
    rows = vars(data),
    cols = vars(chains),
    drop = FALSE,
    scales = "free_y") +
  theme(
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines")
  )

ggsave(
  file = here("figs/overlap_diff_genes.pdf"),
  width = 16,
  height = 9
)


## Original version with counts
# g <- ggplot(
#   all_overlap_df,
#     aes(
#       x = n,
#       fill = var,
#       y = ..count..
#     )
#   ) +
#   geom_bar(position = position_dodge(preserve = "single")) +
#   labs(
#     x = "Number of times gene was identified",
#     y = "Frequency + 1"
#   ) +
#   # scale_y_log10() +
#   theme(
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15)
#   ) +
#   scale_fill_brewer(palette = "Set1", name = "Parameter",
#     limits = c("mu", "delta", "epsilon"),
#     labels = c(bquote(mu[i]), bquote(delta[i]), bquote(epsilon[i]))) +
#   facet_grid(
#     rows = vars(data),
#     cols = vars(chains), scales = "free_y") +
#   theme(legend.position = "bottom")
