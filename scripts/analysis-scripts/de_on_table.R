###############################################################################
##
## Run TestDE
##
###############################################################################

edr <- mclapply(
  seq_len(nrow(df)),
  function(i) {
    cat(i, "/", nrow(df), "\n")
    if (isTRUE(df[[i, "chains"]] == 1)) {
      return(rep(list(NULL), 3))
    }
    ind <- references[["data"]] == df[[i, "data"]]
    chain <- readRDS(df[[i, "file"]])
    if (length(chain) > 1) {
      suppressMessages(
        chain <- Scalability:::combine_subposteriors(
          chain,
          subset_by = "gene",
          mc.cores = 1
        )
      )
    }
    cp <- intersect(c("nu", "s", "phi"), names(chain@parameters))
    chain@parameters[cp] <- lapply(
      chain@parameters[cp],
      function(x) {
        colnames(x) <- gsub("_Batch.*", "", colnames(x))
        x
      }
    )
    chain@parameters <- Scalability:::reorder_params(
      chain@parameters,
      gene_order = rownames(references[[which(ind), "chain"]]),
      cell_order = gsub(
        "_Batch.*",
        "",
        colnames(references[[which(ind), "chain"]])
      )
    )
    nsamples <- nrow(
      references[[which(ind), "chain"]]@parameters[["mu"]]
    )
    chain@parameters <- lapply(
      chain@parameters,
      function(x) {
        x[seq_len(nsamples), ]
      }
    )
    suppressMessages(
      de <- BASiCS_TestDE(
        references[[which(ind), "chain"]],
        chain,
        Plot = FALSE,
        PlotOffset = FALSE
      )
    )
    lapply(
      de@Results,
      function(x) {
        l <- DiffRes(x)[["GeneName"]]
        if (!length(l)) NULL else l
      }
    )
  },
  mc.cores = 2
)


edr_df <- do.call(rbind, edr)
colnames(edr_df) <- c("DiffExp", "DiffDisp", "DiffResDisp")
edr_df <- as.data.frame(edr_df)
df <- cbind(df, edr_df)


df[, c("nDiffExp", "nDiffDisp", "nDiffResDisp")] <- lapply(
  df[, c("DiffExp", "DiffDisp", "DiffResDisp")],
  function(x) sapply(x, length)
)

df[c("pDiffExp", "pDiffDisp", "pDiffResDisp")] <- round(
  df[c("nDiffExp", "nDiffDisp", "nDiffResDisp")] / df[["nGenes"]],
  digits = 3
)




###############################################################################
##
## DE numbers
##
###############################################################################

sdf  <- df[, 
  c("data", "chains", "pDiffExp", "pDiffDisp", "pDiffResDisp")
]
mdf <- melt(sdf, id.vars = c("data", "chains"))
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

mdf <- mdf %>% filter(variable %in% c("mu", "delta"))

advi_mdf <- mdf %>% filter(is.na(chains)) %>% 
  group_by(data, variable, chains) %>% 
  summarize(value = median(value))


ggplot(mdf[!(is.na(mdf$chains) | mdf$chains == 1), ],
  aes(
    x = factor(chains),
    y = value,
    # group = data,
    color = variable
  )
) + 
  geom_quasirandom(
    groupOnX = TRUE,
    dodge.width = 1,
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
    aes(
      yintercept = value,
      color = data,
      linetype = "ADVI")
  ) +
  scale_linetype_manual(name = NULL, labels = "Mean ADVI results", values = 2) +
  facet_wrap(~data, nrow = 2, ncol = 2) +
  scale_x_discrete(name = "Partitions") +
  scale_y_continuous(name = "Portion of genes perturbed", labels = scales::percent) +
  scale_color_brewer(name = "Data", palette = "Set2")

ggsave(here("figs/diffexp_plot.pdf"), width = 12, height = 8)

###############################################################################
##
## DE overlaps
##
###############################################################################

md <- df %>% group_by(chains, data, by)


summ <- function(x) {
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
  mu = summ(DiffExp),
  delta = summ(DiffDisp),
  epsilon = summ(DiffResDisp)
)




overlap_df <- function(mds, i, var) {
  o <- mds[[i, var]]
  if (length(o) > 0 && !all(is.na(names(o)))) {
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
all_overlap_df$chains <- factor(
  paste(all_overlap_df$chains, "chains"),
  levels = c(paste(sort(as.numeric(unique(all_overlap_df$chains))), "chains"), "ADVI")
)
all_overlap_df$chains[is.na(all_overlap_df$chains)] <- "ADVI"
# all_overlap_df <- all_overlap_df[!is.na(all_overlap_df$chains), ]


count_df <- all_overlap_df %>% 
  group_by(n, data, chains, var, nGenes) %>% 
  summarise(count = n())

## Proportions because of differing number of genes
ggplot(
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
# ggplot(
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
