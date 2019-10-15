options(stringsAsFactors = FALSE)
library("devtools")
if (!"package:BASiCS" %in% search()) load_all("../BASiCS")
if (!require("Scalability", quietly = TRUE)) load_all()
library("dplyr")
library("flexmix")
library("transport")
library("ggplot2")
library("ggrepel")
library("RColorBrewer")
library("ComplexHeatmap")
theme_set(theme_bw())



dir.create("figs", showWarnings = FALSE)
data <- list.files("datasets", pattern = "Data")

metadata <- data.frame(
  file = data,
  dataset = gsub("(.*)Data.*", "\\1", data),
  chains = gsub(".*Data_by_(cell|gene)_c_?(\\d+).*", "\\2", data),
  by = gsub(".*Data_by_(cell|gene)_c_?(\\d+).*", "\\1", data),
  seed = gsub(".*Data_by_(cell|gene)_c_?\\d+_s_?(\\d+).*", "\\2", data),
  method = gsub(".*Data_by_(cell|gene)_c_?\\d+_s_?\\d+_?(.*)\\.rds", "\\2", data)
)

ind_keep <- na.omit(union(which(metadata$chains == "01")[1], which(metadata$chains != "01")))
metadata <- metadata[ind_keep, ]
metadata <- metadata[metadata$dataset == "zeisel", ]
metadata$method[metadata$method == ""] <- "none"
metadata$by[metadata$method == "none"] <- "none"
# metadata <- metadata[metadata$method %in% c("none", "pien"), ]
metadata$chains <- as.numeric(as.character(metadata$chains))

## For doing DE plot
# metadata <- metadata[metadata$by != "gene", ]

c1 <- readRDS(file.path("datasets", metadata$file[[1]]))


# metadata$objects <- lapply(seq_len(nrow(metadata)), function(i) {
#   print(i)
#   readRDS(file.path("datasets", metadata[i, "file"]))
# }
# # , mc.cores = 16
# )




gene_order <- colnames(c1@parameters[["mu"]])
cell_order <- colnames(c1@parameters[["s"]])


# for (i in seq_along(metadata$objects)) {
#   x <- metadata$objects[[i]]
#   metadata$objects[[i]] <- reorder_chain(x, gene_order, cell_order)
# }


# i <- which(metadata$chains == 1)
# c1 <- metadata$objects[[i]]
# for (j in setdiff(seq_along(metadata$objects), i)) {
#   metadata$objects[[j]] <- offset_correct(c1, metadata$objects[[j]])
# }


diffexp <- function(foo) {
  !foo %in% c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")
}  

de_res <- mclapply(seq_len(nrow(metadata)),
  function(j) {
    if (i == j) {
      return(NULL)
    }
    print(j)
    c2 <- readRDS(file.path("datasets", metadata[j, "file"]))
    c2 <- reorder_chain(c2, gene_order, cell_order)
    # c2 <- metadata$objects[[j]]
    de <- BASiCS_TestDE(
      c1,
      c2,
      Plot = FALSE,
      PlotOffset = FALSE
    )
    data.frame(
      Ind = j,
      "DiffMean" = sum(diffexp(de$TableMean$ResultDiffMean)),
      "DiffDisp" = sum(diffexp(de$TableDisp$ResultDiffDisp)),
      "DiffResDisp" = sum(diffexp(de$TableResDisp$ResultDiffResDisp))
    )
  },  mc.cores = 16)
de_res <- de_res[-1]
de_res <- do.call(rbind, de_res)
de_res <- cbind(de_res, metadata[-1, setdiff(colnames(metadata), c("objects", "summaries"))])


mdf <- reshape2::melt(de_res, measure.var = c("DiffMean", "DiffDisp", "DiffResDisp"))
mdf$variable <- gsub("DiffMean", "mu", mdf$variable)
mdf$variable <- gsub("DiffDisp", "delta", mdf$variable)
mdf$variable <- gsub("DiffResDisp", "epsilon", mdf$variable)
mdf$variable <- factor(mdf$variable)


mdf$method <- gsub("var$", " (variance weighting)", mdf$method)
mdf$method <- gsub("n$", " (n weighting)", mdf$method)
mdf$method <- sub("([[:alpha:]])", "\\U\\1", mdf$method, perl = TRUE)
mdf$method <- gsub("Pie", "PIE", mdf$method, perl = TRUE)


mean_diff <- mdf %>% 
  dplyr::filter(by == "cell") %>%
  dplyr::group_by(variable, by, chains, method) %>%
  dplyr::summarise(median = median(value), min = min(value), max = max(value))
mean_diff$variable <- factor(mean_diff$variable, levels = c("mu", "delta", "epsilon"))

facet_names <- c(
  mu = "mu[i]",
  delta = "delta[i]",
  epsilon = "epsilon[i]"
)
mean_diff$method <- gsub("(", "\n(", mean_diff$method, fixed = TRUE)

g <- ggplot(
    mean_diff, 
    aes(x = factor(chains), ymin = min, y = median, ymax = max, colour = method)
  ) + 
  geom_pointrange(fatten = 1.5, position = position_dodge(width = 0.75)) +
  labs(x = "Number of partitions", y = "Median number of genes") +
  scale_colour_brewer(palette = "Set2", name = "Method") +
  facet_wrap(~variable, labeller = as_labeller(facet_names, label_parsed)) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2,byrow = TRUE))

ggsave(g, file = "figs/de_by_methods.pdf", width = 5, height = 4)



min_median <- mean_diff %>% 
  dplyr::group_by(variable, by, chains) %>%
  dplyr::filter(median == min(median))

tab <- table(min_median$method)



de_res_f <- de_res[de_res$method %in% c("pien", "none"), ]
metadata_f <- metadata[metadata$method %in% c("pien", "none"), ]

mdf <- reshape2::melt(de_res_f, measure.var = c("DiffMean", "DiffDisp", "DiffResDisp"))
mdf$variable <- gsub("DiffMean", "mu", mdf$variable)
mdf$variable <- gsub("DiffDisp", "delta", mdf$variable)
mdf$variable <- gsub("DiffResDisp", "epsilon", mdf$variable)
mdf$variable <- factor(mdf$variable, levels = c("mu", "delta", "epsilon"))





################################################################################
#
# Differential expression figure
#
################################################################################
mdf$by <- paste("Partitioned by", mdf$by)

mdf[mdf$variable == "mu", "variable"] <- bquote(mu[i])
mdf[mdf$variable == "delta", "variable"] <- bquote(delta[i])
mdf[mdf$variable == "epsilon", "variable"] <- bquote(epsilon[i])


g <- ggplot(
    mdf,
    aes(x = chains, y = value + 1, colour = variable
      # , shape = by
      )
    ) + 
  geom_jitter(width = 0.2, height = 0, size = 1) +
  geom_segment(
    data = data.frame(
      x = c(1, 0.4, 1),
      xend = c(200, 200, 200),
      y = c(1.67, 1.67, 71.67) + 1,
      yend = c(1.67, 1.67, 71.67) + 1,
      color = c("mu", "delta", "epsilon"),
      size = c("mean ADVI results", "mean ADVI results", "mean ADVI results")
    ),
    linetype = "dashed", 
    mapping = aes(
      x = x,
      xend = xend,
      y = y, 
      yend = yend,
      color = color,
      size = size
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(trans = "log2", breaks = sort(unique(metadata_f$chains))) + 
  scale_y_log10() +
  scale_size_manual(name = NULL, limits = "mean ADVI results", values = 1) +
  coord_cartesian(xlim = c(1.7, 150)) +
  guides(linetype = FALSE) +
  scale_colour_brewer(palette = "Set1", name = "Parameter",
    limits = c("mu", "delta", "epsilon"),
    labels = c(bquote(mu[i]), bquote(delta[i]), bquote(epsilon[i]))) +
  scale_linetype_discrete(name = NULL) +
  facet_wrap(~by) +
  labs(x = "Number of partitions", y = "Number of detected genes + 1") +
  theme(legend.position = "bottom")

ggsave(g, file = "figs/de_vs_chains.pdf", width = 6, height = 4)



fpr_table <- dplyr::filter(de_res_f) %>% 
  dplyr::group_by(chains, by) %>% 
  dplyr::summarise(
    DiffMean = mean(DiffMean), 
    DiffDisp = mean(DiffDisp), 
    DiffResDisp = mean(DiffResDisp)
  )
fpr_table <- fpr_table[order(fpr_table$by), ]
colnames(fpr_table) <- c(
  "# Chains", 
  "Partitioned by",
  "Expression", 
  "Over-dispersion", 
  "Residual over-dispersion")
fpr_table[, 1] <- lapply(fpr_table[, 1], as.integer)


print(xtable::xtable(fpr_table, 
  digits = 1,
  label = "table::DiffResDC",
  caption = "Mean number of genes highlighted as having differential expression
    and differential over-dispersion genes in different runs of divide and 
    conquer MCMC, relative to MCMC using the full data.",
  ), 
  include.rownames = FALSE, 
  file = "tables/divide_and_conquer_TestDE.tex"
)



################################################################################
#
# Trends
#
################################################################################



get_file <- function(by, nchains, seed) {
  paste0(
    "/exports/eddie/scratch/s1372510/Benchmarks/zeiselData/", 
    by, "/", 
    sprintf("%02d", nchains),
    "_chains_seed_", 
    seed, 
    ".rds")
}

g1 <- plot_fits(readRDS(get_file("gene", 2, 42)), c1, alpha = 0.6)
g2 <- plot_fits(readRDS(get_file("gene", 128, 42)), c1, alpha = 0.6)

ggsave(g1 + labs(x = bquote(log(mu[i])), y = bquote(log(delta[i]))),
  file = "figs/trend_2.pdf", 
  width = 4, height = 3)

ggsave(g2 + labs(x = bquote(log(mu[i])), y = bquote(log(delta[i]))),
  file = "figs/trend_128.pdf", 
  width = 4, height = 3)






################################################################################
#
# MA plots and some violin
#
################################################################################

s1 <- Summary(c1)
for (param in c("mu", "delta", "epsilon")) {
  i2 <- which(metadata_f$by == "cell" & metadata_f$chains == 32)[[1]]
  c2 <- readRDS(file.path("datasets", metadata_f[i2, "file"]))
  c2 <- reorder_chain(c2, gene_order, cell_order)
  c2 <- offset_correct(c1, c2)
  s2 <- Summary(c2)

  if (param == "epsilon") {
    yl <- bquote(epsilon[i]^F - epsilon[i]^DC)
  } else if (param == "delta") {
    yl <- bquote(log2(delta[i]^F / delta[i]^DC))
  } else if (param == "mu") {
    yl <- bquote(log2(mu[i]^F / mu[i]^DC))
  }
  g <- plot_param_ma(
    s1, 
    "No partitioning", 
    s2, 
    paste0(metadata_f[i2, "chains"], " partitions (", metadata_f[i2, "by"], ")"),
    param = param, log = param != "epsilon"
  ) + 
    ggtitle(NULL) + 
    labs(y = yl, x = bquote(mu[i]))
  desc <- paste(metadata_f[i2, "by"], metadata_f[i2, "chains"], sep = "_")
  ggsave(
    g,
    file = paste0("figs/cell_32_ma_", param, ".pdf"),
    width = 3,
    height = 3
  )
  i2 <- which(metadata_f$by == "gene" & metadata_f$chains == 128)[[1]]
  c2 <- readRDS(file.path("datasets", metadata_f[i2, "file"]))
  c2 <- reorder_chain(c2, gene_order, cell_order)
  c2 <- offset_correct(c1, c2)
  s2 <- Summary(c2)
  g <- plot_param_ma(
    s1, 
    "No partitioning", 
    s2, 
    paste0(metadata_f[i2, "chains"], " partitions (", metadata_f[i2, "by"], ")"),
    param = param, log = param != "epsilon"
  ) + ggtitle(NULL) + ylab(yl) + xlab(bquote(mu[i]))
  desc <- paste(metadata_f[i2, "by"], metadata_f[i2, "chains"], sep = "_")
  ggsave(
    g,
    file = paste0("figs/gene_128_ma_", param, ".pdf"),
    width = 3,
    height = 3
  )
}


################################################################################
#
# computational time
#
################################################################################


times <- readRDS("datasets/times.rds")
dir <- "/exports/eddie/scratch/s1372510/Benchmarks/"
time_metadata <- data.frame(
  file = names(times),
  data = gsub(
    paste0(dir, "/(\\w+)Data/.*"), 
    "\\1", 
    names(times)),
  by = gsub(
    paste0(dir, "/\\w+/(\\w+)/.*"), 
    "\\1", 
    names(times)),
  chains = gsub(
    paste0(dir, "/\\w+/\\w+/(\\d+).*"), 
    "\\1", 
    names(times)),
  seed = gsub(
    paste0(dir, "/\\w+/\\w+/\\d+_chains_seed_(\\d+)_time\\.rds"), 
    "\\1", 
    names(times)),
  time = sapply(times, function(x) x[[3]])
)

time_metadata$chains <- as.numeric(as.character(time_metadata$chains))
# n <- unique(time_metadata$chains)
# time_metadata$chains <- factor(time_metadata$chains, levels = n[order(as.numeric(n))])
time_metadata$data <- sub("([[:alpha:]])", "\\U\\1", time_metadata$data, perl = TRUE)
time_metadata$by <- sub("([[:alpha:]])", "\\U\\1", time_metadata$by, perl = TRUE)


tmax <- max(time_metadata$time)

averages <- time_metadata %>% group_by(by, chains) %>% summarise(time = mean(time))
averages <- as.data.frame(averages)

## sub in 32, 64, 128 times because of 16 core limitation
averages[averages$chains == 32, "time"] <- mean(readRDS("datasets/t32.rds"))
averages[averages$chains == 64, "time"] <- mean(readRDS("datasets/t64.rds"))
averages[averages$chains == 128, "time"] <- mean(readRDS("datasets/t128.rds"))


averages$n <- NA
averages$n[averages$by == "Gene"] <- ncol(c1@parameters[["mu"]]) / 
  averages[averages$by == "Gene", "chains", drop = TRUE]
averages$n[averages$by == "Cell"] <- ncol(c1@parameters[["s"]]) / 
  averages[averages$by == "Cell", "chains", drop = TRUE]

out <- averages %>% 
  group_by(by) %>% 
  dplyr::mutate(n_str = paste0(round(n), " ", tolower(by), "s"))

out_f <- dplyr::filter(out, by == "Gene")
out <- as.data.frame(out)
out_f <- as.data.frame(out_f)

out_f[1:6, "n_str"] <- paste(out[1:6, "n_str"], out_f[1:6, "n_str"], sep = ",\n")


fs <- 2.25

g <- ggplot(out_f, 
    aes(
      x = chains, 
      y = time / 3600, 
    )
  ) +
  geom_point() +
  geom_hline(yintercept = 12255 / 3600, lty = "dashed", colour = "grey60") + 
  annotate(x = 1, y = 12255 / 3600 * 1.2, size = fs,
    label = "Time taken for ADVI", geom = "text", hjust = 0) +
  geom_hline(yintercept = tmax / (10 * 3600), lty = "dashed", colour = "grey60") + 
  annotate(x = 1, y = tmax / (10 * 3600) * 1.2, size = fs,
    label = "10x speedup", geom = "text", hjust = 0) +
  geom_hline(yintercept = tmax / (100 * 3600), lty = "dashed", colour = "grey60") + 
  annotate(x = 1, y = tmax / (100 * 3600) * 1.2, size = fs,
    label = "100x speedup", geom = "text", hjust = 0) +
  geom_line() +
  geom_text(
    data = out_f, 
    hjust = 0.4,
    size = fs,
    mapping = aes(
      x = chains,
      y = 25,
      label = n_str),
    show.legend = FALSE
  ) +
  labs(x = "Number of partitions", y = "Time") +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64, 128)) +
  scale_y_continuous(
    trans = "log10", 
    breaks = c(0.16666, 0.5, 1, 2, 5, 15), 
    labels = c("10 min", "30 min", "1 hr", "2 hr", "5 hr", "15 hr")
  )

ggsave(g, file = "figs/time_plot.pdf", width = 6, height = 4)






################################################################################
#
# Plot minimum ESS
#
################################################################################





calc_median_ess <- function(i) {
  print(i)
  chain <- readRDS(file.path("datasets", metadata[i, "file"]))
  data.frame(
    mu = median(BASiCS_effectiveSize(chain, "mu"), na.rm = TRUE),
    delta = median(BASiCS_effectiveSize(chain, "delta"), na.rm = TRUE),
    epsilon = median(BASiCS_effectiveSize(chain, "epsilon"), na.rm = TRUE)
  )
}

ess <- mclapply(
  seq_len(nrow(metadata)), 
  calc_median_ess, 
  mc.cores = 16
)

essdf <- do.call(rbind, ess)
essdf <- cbind(essdf, metadata[, c("chains", "by", "method")])


essmdf <- reshape2::melt(essdf, measure.vars = c("mu", "delta", "epsilon"))
essmdf <- essmdf[essmdf$by != "none", ]
essmdf$by <- paste("Partitioned by", essmdf$by)
essmdf$method <- gsub("pie", "PIE", essmdf$method)
essmdf$method <- gsub("consensus", "Consensus", essmdf$method)
essmdf$method <- gsub("n$", " (n weighting)", essmdf$method)
essmdf$method <- gsub("var$", " (variance weighting)", essmdf$method)
essmdf$method <- gsub("(", "\n(", essmdf$method, fixed = TRUE)

g <- ggplot(
    essmdf, 
    aes(x = chains, y = value, colour = variable)) + 
  geom_hline(yintercept = 1000, lty = "dashed", colour = "grey60") +
  geom_point() +
  scale_x_continuous(trans = "log2", breaks = c(2, 4, 8, 16, 32, 64, 128)) +
  scale_colour_brewer(palette = "Set1", name = "Parameter",
    limits = c("mu", "delta", "epsilon"),
    labels = c(bquote(mu[i]), bquote(delta[i]), bquote(epsilon[i]))) +
  labs(x = "Number of partitions", y = "Median effective sample size") +
  facet_grid(rows = vars(by), cols = vars(method)) + 
  theme(legend.position = "bottom")

ggsave(g, file = "figs/median_ess.pdf", width = 7, height = 4)


pie_ess_chain <- readRDS(
  file.path(
    "datasets",
    metadata_f[which(metadata_f$method=="pien" & metadata_f$by == "cell")[[1]], "file"]
  )
)
pie_ess_chain2 <- pie_ess_chain
pie_ess_chain2@parameters <- lapply(pie_ess_chain2@parameters, 
  function(x) {
    apply(x, 2, function(col) sample(col, length(col)))
  })

g1 <- BASiCS_diagHist(pie_ess_chain)
g2 <- BASiCS_diagHist(pie_ess_chain2)

ggsave(
  g1 +
    scale_x_log10() +
    scale_y_log10(),
  file = "figs/pie_ess.pdf", width = 4, height = 4
)
ggsave(
  g2 +
    scale_x_log10() +
    scale_y_log10(),
  file = "figs/pie_ess_fixed.pdf", width = 4, height = 4
)



do_ess_hist <- function(i) {
  chain <- readRDS(file.path("datasets", metadata_f[i, "file"]))
  BASiCS_diagHist(chain)
}

do_ess_hist(which.min(essdf$mu))
do_ess_hist(which.min(essdf$delta))
do_ess_hist(which.min(essdf$epsilon))


################################################################################
#
# density plot to demonstrate technique
#
################################################################################
do_density_plot <- function(param, by, nchains, seed, ngenes = 5, log = TRUE) {
  ind <- which(
    metadata_f$by == by & 
      metadata_f$chains == nchains & 
      metadata_f$seed == seed
  )
  dc <- readRDS(file.path("datasets", metadata_f$file[[ind]]))
  dc <- reorder_chain(dc, gene_order, cell_order)
  cd <- offset_correct(c1, dc)
  dcc <- readRDS(get_file(by, nchains, seed))
  if (by == "cell") {
    dcc[] <- lapply(dcc, function(x) offset_correct(c1, x))
  }
  discordant <- abs(colMedians(dc@parameters[[param]]) - colMedians(c1@parameters[[param]]))
  names(discordant) <- colnames(dc@parameters[[param]])
  discordant <- sort(discordant, decreasing = TRUE)[1:ngenes]
  most_discordant <- names(discordant)
  parameter_densities(
    param,
    reference_chain = c1,
    chains = dcc,
    collapsed_chain = dc,
    columns = most_discordant,
    log = log
  )
}



g1 <- do_density_plot("mu", "cell", 2, 7, 5)
ggsave(
  g1 +
    ggtitle(NULL) +
    xlab("Value"),
  file = "figs/cell_densities.pdf", width = 4, height = 6
)
g2 <- do_density_plot("mu", "gene", 2, 7, 4)
ggsave(
  g2 +
    ggtitle(NULL) +
    xlab("Value"),
  file = "figs/gene_densities.pdf", width = 4, height = 6
)
# gd <- do_density_plot("delta", "cell", 8, 7, 5)
# g3 <- do_density_plot("delta", "gene", 8, 14, 5)



################################################################################
#
# identity/overlap of genes
#
################################################################################

de_objs <- mclapply(seq_len(nrow(metadata_f)), function(j) {
  if (i == j) {
    return(NULL)
  }
  print(j)
  c2 <- readRDS(file.path("datasets", metadata_f[j, "file"]))
  c2 <- reorder_chain(c2, gene_order, cell_order)
  c2 <- offset_correct(c1, c2)
  # c2 <- metadata_f$objects[[j]]
  BASiCS_TestDE(
    c1, 
    c2, 
    Plot = FALSE, 
    PlotOffset = FALSE
  )
}, mc.cores = 16)

get_genes <- function(i, which) {
  table <- de_objs[[i]][[paste0("Table", which)]]
  ind <- which(diffexp(table[[paste0("ResultDiff", which)]]))
  table[ind, "GeneName"]
}

lm <- mclapply(
  seq_len(nrow(metadata_f)),
  function(i) get_genes(i, "Mean"),
  mc.cores = 16
)
ld <- mclapply(
  seq_len(nrow(metadata_f)),
  function(i) get_genes(i, "Disp"),
  mc.cores = 16
)
le <- mclapply(
  seq_len(nrow(metadata_f)),
  function(i) get_genes(i, "ResDisp"),
  mc.cores = 16
)
names(lm) <- seq_along(lm)
names(ld) <- seq_along(ld)
names(le) <- seq_along(le)


metadata_f$mu_genes <- lm
metadata_f$delta_genes <- ld
metadata_f$epsilon_genes <- le

md <- metadata_f %>% group_by(by, chains)


summ <- function(x) {
  list(
    sapply(Reduce(union, x), 
      function(y) {
        sum(sapply(x, function(z) y %in% z))
      }
    )
  )
}
mds <- md %>% dplyr::summarise(
  mu = summ(mu_genes),
  delta = summ(delta_genes),
  epsilon = summ(epsilon_genes)
)

get_overlaps <- function(mds, split, chain, column) { 
  n <- mds[mds$by == split & mds$chains == chain, column]
  n[[1]][[1]]
}


overlap_df <- function(mds, by, chains, var) {
  o <- get_overlaps(mds, by, chains, var)
  if (length(o)) {
    data.frame(
      by = by,
      chains = chains,
      var = var,
      n = o
    )      
  } else {
    as.data.frame(
      matrix(
        nrow = 0, 
        ncol = 4, 
        dimnames = list(NULL, c("by", "chains", "var", "n"))
      )
    )
  }
}

all_overlap_df <- lapply(seq_len(nrow(mds)), function(i) {
  d1 <- overlap_df(mds, mds[[i, "by"]], mds[[i, "chains"]], "mu")
  d2 <- overlap_df(mds, mds[[i, "by"]], mds[[i, "chains"]], "delta")
  d3 <- overlap_df(mds, mds[[i, "by"]], mds[[i, "chains"]], "epsilon")
  rbind(d1, d2, d3)
})
all_overlap_df <- do.call(rbind, all_overlap_df)
all_overlap_df$n <- factor(all_overlap_df$n)
all_overlap_df$chains <- factor(
  paste(all_overlap_df$chains, "chains"),
  levels = paste(sort(unique(all_overlap_df$chains)), "chains")
)
# all_overlap_df$var <- factor(all_overlap_df$var, levels = c("mu", "delta", "epsilon"))
all_overlap_df$by <- paste("Partitioned by", all_overlap_df$by)
g <- ggplot(all_overlap_df, aes(x = n, fill = var, y = ..count.. + 1)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  labs(
    x = "Number of times gene was identified",
    y = "Frequency + 1"
  ) +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  scale_fill_brewer(palette = "Set1", name = "Parameter",
    limits = c("mu", "delta", "epsilon"),
    labels = c(bquote(mu[i]), bquote(delta[i]), bquote(epsilon[i]))) +
  facet_grid(rows = vars(by), cols = vars(chains), scales = "free_y") +
  theme(legend.position = "bottom")


ggsave(
  g, 
  file = "figs/overlap_diff_genes.pdf",
  width = 10,
  height = 6)



################################################################################
#
# VB stuff
#
################################################################################

vb_objects <- readRDS("outputs/vb_objects.rds")
elbos <- lapply(vb_objects, function(x) x[[1]])
vb_objects <- lapply(vb_objects, function(x) x[[2]])


## Convergence first
parse_elbo <- function(c) {
  c <- gsub("Chain 1:\\s+", "", c)
  elbo <- c[(grep("Begin stochastic", c) + 1):(grep("Drawing", c) - 2)]
  elbo[-c(1, grep("CONVERGED", elbo))] <- paste(
    elbo[-c(1, grep("CONVERGED", elbo))],
    "NOTCONVERGED")
  elbo <- gsub("(MEDIAN )?ELBO CONVERGED", "CONVERGED", elbo)
  elbo <- strsplit(elbo, "\\s+")
  elbo <- do.call(rbind, elbo)
  colnames(elbo) <- elbo[1, ]
  elbo <- elbo[-1, ]
  elbo <- as.data.frame(elbo, stringsAsFactors=FALSE)
  elbo[, 1:4] <- lapply(elbo[, 1:4], as.numeric)
  elbo
}

elbo_dfs <- lapply(elbos, parse_elbo)
elbo_dfs <- lapply(seq_along(elbo_dfs), function(i) {
  x <- elbo_dfs[[i]]
  x$chain <- i
  x
})
elbo_df <- do.call(rbind, elbo_dfs)


g <- ggplot(elbo_df, aes(x = iter, y = ELBO, colour = factor(chain))) + 
  geom_line(alpha = 0.8) +
  labs(x = "Iteration") +
  scale_colour_brewer(palette = "Dark2") +
  guides(colour = FALSE)

ggsave(g, file = "figs/ELBO_convergence.pdf", width = 5, height = 3.5)



vb_b_objects <- lapply(vb_objects, function(x) {
  xe <- extract(x)
  parameters <- list(
    mu = xe$mu,
    delta = xe$delta,
    epsilon = xe$epsilon,
    s = xe$s,
    nu = xe$nu,
    beta = xe$beta,
    theta = as.matrix(xe$theta),
    phi = xe$phi
  )
  gp <- c("mu", "delta", "epsilon")
  cp <- c("s", "nu", "phi")
  parameters[gp] <- lapply(gp, function(x) {
    colnames(parameters[[x]]) <- colnames(c1@parameters[["mu"]])
    parameters[[x]]
  })
  parameters[cp] <- lapply(cp, function(x) {
    colnames(parameters[[x]]) <- colnames(c1@parameters[["nu"]])
    parameters[[x]]
  })

  new("BASiCS_Chain", parameters = parameters)
})


vb_res <- lapply(seq_along(vb_b_objects), function(i) {
  vb <- vb_b_objects[[i]]
  de <- BASiCS_TestDE(
    vb, 
    c1, 
    GroupLabel1 = "VB", 
    GroupLabel2 = "MCMC", 
    Plot = FALSE,
    PlotOffset = FALSE
  )
  data.frame(
    Ind = i,
    "DiffMean" = sum(diffexp(de$TableMean$ResultDiffMean)),
    "DiffDisp" = sum(diffexp(de$TableMean$ResultDiffMean)),
    "DiffResDisp" = sum(diffexp(de$TableResDisp$ResultDiffResDisp))
  )
})

vb_resdf <- do.call(rbind, vb_res)

colnames(vb_resdf) <- c(
  "Instance", 
  "Expression", 
  "Over-dispersion", 
  "Residual over-dispersion"
)
vb_resdf <- rbind(vb_resdf, c("Mean", round(colMeans(vb_resdf[, 2:4]), digits = 2)))

print(
  xtable::xtable(
    vb_resdf, 
    label = "table::DiffResVB",
    caption = paste(
      "Number of genes highlighted as having differential expression and", 
      "differential (residual) over-dispersion in different runs of ADVI,
      relative to MCMC."
    )
  ), 
  include.rownames = FALSE, 
  file = "tables/vb_TestDE.tex"
)



hpd_width <- function(summary, param) {
  abs(summary@parameters[[param]][, 2] - summary@parameters[[param]][, 3])
}


hpd_violin <- function(b, vb_b_objects, param) {
  df <- data.frame(MCMC = extract_hpd_interval(Summary(b), param))
  df <- cbind(df, 
    do.call(
      cbind, 
      lapply(vb_b_objects, function(x) extract_hpd_interval(Summary(x), param))
    )
  )
  colnames(df)[-1] <- paste("ADVI run", seq_len(ncol(df) - 1))
  ggplot(melt(as.matrix(df)), aes(x = Var2, y = value)) + 
    geom_violin(fill = "grey80") + 
    geom_boxplot(width = 0.15, fill = "grey80", outlier.colour = NA) +
    labs(x = NULL) +
    scale_y_log10(name = paste0("HPD width (", param, ")"))
}
g1 <- hpd_violin(c1, vb_b_objects, "mu")
g2 <- hpd_violin(c1, vb_b_objects, "delta")
g3 <- hpd_violin(c1, vb_b_objects, "epsilon")

ggsave(g1, file = "figs/hpd_violin_mu.pdf", width = 6, height = 2)
ggsave(g2, file = "figs/hpd_violin_delta.pdf", width = 6, height = 2)
ggsave(g3, file = "figs/hpd_violin_epsilon.pdf", width = 6, height = 2)



hpd_stratified_violin <- function(b, vb, param) {
  s1 <- Summary(b)
  df <- data.frame(MCMC = extract_hpd_interval(s1, param))
  df <- cbind(df, 
    ADVI = extract_hpd_interval(Summary(vb_b_objects[[1]]), param),
    Mu = s1@parameters[["mu"]][, 1]
  )
  if (param == "mu") {
    df[c("MCMC", "ADVI")] <- df[c("MCMC", "ADVI")] / df[["Mu"]]
    yl <- bquote("Normalised HPD width" ~(mu[i]))
  } else if (param == "delta") {
    yl <- bquote("HPD width" ~(delta[i]))
  } else {
    yl <- bquote("HPD width" ~(epsilon[i]))
  }
  df$Quantile <- cut(
    df$Mu, 
    quantile(df$Mu, probs = seq(0, 1, length.out = 11)),
    include.lowest = TRUE)
  df$Quantile <- as.numeric(df$Quantile)
  df$Mu <- NULL
  mdf_v <- reshape2::melt(df, id.var = "Quantile")

  ggplot(mdf_v, 
    aes(
      x = factor(Quantile), 
      y = value, 
      colour = variable, 
      fill = variable)) + 
    geom_violin() + 
    scale_y_log10() +
    labs(x = bquote(Decile ~(mu[i])), y = yl)
}



g1 <- hpd_stratified_violin(c1, vb_b_objects[[1]], "mu")
g2 <- hpd_stratified_violin(c1, vb_b_objects[[1]], "delta")
g3 <- hpd_stratified_violin(c1, vb_b_objects[[1]], "epsilon")


ggsave(g1, file = "figs/stratified_hpd_violin_mu.pdf", width = 6, height = 4)
ggsave(g2, file = "figs/stratified_hpd_violin_delta.pdf", width = 6, height = 4)
ggsave(g3, file = "figs/stratified_hpd_violin_epsilon.pdf", width = 6, height = 4)

################################################################################
#
# wasserstein/Kullback-leibler
#
################################################################################


## distance/divergence for a chain param by param
metadata_f$wasserstein <- vapply(
  seq_len(nrow(metadata_f)),
  function(i) {

    c2 <- readRDS(file.path("datasets", metadata_f[i, "file"]))
    c2 <- reorder_chain(c2, gene_order, cell_order)
    c2 <- offset_correct(c1, c2)
    wd <- vapply(
      colnames(c1@parameters[[param]]), 
      function(col) {
        wasserstein1d(
          c2@parameters[[param]][, col], 
          c1@parameters[[param]][, col]
        )
      }, numeric(1)
    )
    mean(wd)
  }, 
  numeric(1)
)

extract_dist <- function(param, fun) {
  vapply(
    seq_len(nrow(metadata_f)),
    function(i) {
      kl <- vapply(
        c2 <- readRDS(file.path("datasets", metadata_f[i, "file"]))
        c2 <- reorder_chain(c2, gene_order, cell_order)
        c2 <- offset_correct(c1, c2)
        colnames(c1@parameters[[param]]), 
        function(col) {
          fun(
            c2@parameters[[param]][, col],
            c1@parameters[[param]][, col])
        }, numeric(1)
      )
      mean(kl)
    }, numeric(1)
  )
}

getKL <- function(a, b) {
  KLdiv(
    cbind(
      col = a, 
      ref = b
    )
  )["col", "ref"]
}



metadata_f$kullback_leibler_mu <- extract_dist("mu", getKL)
metadata_f$kullback_leibler_delta <- extract_dist("delta", getKL)
metadata_f$kullback_leibler_epsilon <- extract_dist("epsilon", getKL)

metadata_f$wasserstein_mu <- extract_dist("mu", wasserstein1d)
metadata_f$wasserstein_delta <- extract_dist("delta", wasserstein1d)
metadata_f$wasserstein_epsilon <- extract_dist("epsilon", wasserstein1d)


ggplot(metadata_f[metadata_f$by != "none", ], aes(x = chains)) + 
  geom_point(aes(y = kullback_leibler_mu, color = "mu")) + 
  geom_point(aes(y = kullback_leibler_delta, color = "delta")) + 
  # geom_point(aes(y = kullback_leibler_epsilon, color = "epsilon")) + 
  labs(x = "Number of partitions", y = "Mean KL divergence (1D)") +
  scale_x_continuous(trans = "log2") +
  scale_colour_manual(
    values = brewer.pal(3, "Set1"), 
    limits = c("mu", "delta", "epsilon")
  ) +
  facet_wrap(~by)

ggplot(metadata_f[metadata_f$by != "none", ], aes(x = chains)) + 
  geom_point(aes(y = wasserstein_mu, color = "mu")) + 
  geom_point(aes(y = wasserstein_delta, color = "delta")) + 
  geom_point(aes(y = wasserstein_epsilon, color = "epsilon")) + 
  labs(x = "Number of partitions", y = "Mean Wasserstein distance (1D)") +
  scale_x_continuous(trans = "log2") +
  scale_colour_manual(
    values = brewer.pal(3, "Set1"), 
    limits = c("mu", "delta", "epsilon")
  ) +
  facet_wrap(~by)
