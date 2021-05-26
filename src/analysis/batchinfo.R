files_b <- list.files(
  list.files(
    list.files("outputs/batchinfo", full.names = TRUE),
    full.names = TRUE
    ),
  full.names = TRUE
)

chains_b <- lapply(files_b, readRDS)
desc_b <- data.frame(
  data = gsub("outputs/\\w+/\\d/(\\w+)/.*$", "\\1", files_b),
  type = gsub("outputs/\\w+/\\d/\\w+/(.*).rds$", "\\1", files_b)
)


bde1 <- BASiCS_TestDE(
  chains_b[[1]],
  chains_b[[2]],
  GroupLabel1 = "Batch",
  GroupLabel2 = "No batch",
  EFDR_M = NULL,
  EFDR_D = NULL,
  EFDR_R = NULL
)

bde2 <- BASiCS_TestDE(
  chains_b[[3]],
  chains_b[[4]],
  GroupLabel1 = "Batch",
  GroupLabel2 = "No batch",
  EFDR_M = NULL,
  EFDR_D = NULL,
  EFDR_R = NULL
)
g1 <- BASiCS_PlotDE(bde1, Plots = c("MAPlot", "VolcanoPlot"))
g2 <- BASiCS_PlotDE(bde2, Plots = c("MAPlot", "VolcanoPlot"))
ggsave(g1, file = "figs/de_batch_tung.pdf", width = 12, height = 10)
ggsave(g2, file = "figs/de_batch_zeisel.pdf", width = 12, height = 10)


ess_b <- lapply(c("mu", "delta", "epsilon"), 
  function(param) {
    lapply(seq_along(chains_b),
      function(i) {
        x <- chains_b[[i]]
        cbind(
          ess = effectiveSize(x@parameters[[param]]),
          gene = colnames(x@parameters[[param]]),
          param = param,
          desc_b[i, ]
        )
      }
    )
  }
)

ess_dfs <- lapply(ess_b, function(x) do.call(rbind, x))
ess_df_b <- do.call(rbind, ess_dfs)
ess_df_b$data <- gsub("tung", "Tung", ess_df_b$data)
ess_df_b$data <- gsub("zeisel", "Zeisel", ess_df_b$data)
ess_df_b[["param"]] <- factor(ess_df_b[["param"]],
  levels = c("mu", "delta", "epsilon")
)

ggplot(ess_df_b, aes(x = param, y = ess, color = type, fill = type)) + 
  geom_violin(alpha = 0.5) +
  facet_wrap(~data) +
  scale_y_log10(name = "Effective sample size") +
  scale_x_discrete(
    name = "Parameter",
    breaks = c("mu", "delta", "epsilon")
  ) +
  scale_fill_brewer(
    name = "Model",
    aesthetics = c("fill", "color"),
    breaks = c("batch", "nobatch"), 
    labels = c("Using batch", "Ignoring batch"),
    palette = "Dark2"
  )

ggsave("figs/ess_batch.pdf", width = 6, height = 4)
