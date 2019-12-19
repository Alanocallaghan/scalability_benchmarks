files_id <- list.files(
  list.files(
    list.files("outputs/identifiable", full.names = TRUE),
    full.names = TRUE
    ),
  full.names = TRUE
)

chains_id <- lapply(files_id, readRDS)
desc_id <- data.frame(
  data = gsub("outputs/\\w+/\\d/(\\w+)/.*$", "\\1", files_id),
  type = gsub("outputs/\\w+/\\d/\\w+/(.*).rds$", "\\1", files_id)
)


ide1 <- BASiCS_TestDE(
  chains_id[[1]],
  chains_id[[2]],
  GroupLabel1 = "ID",
  GroupLabel2 = "Non-ID",
  EFDR_M = NULL,
  EFDR_D = NULL,
  EFDR_R = NULL
)

ide2 <- BASiCS_TestDE(
  chains_id[[3]],
  chains_id[[4]],
  GroupLabel1 = "ID",
  GroupLabel2 = "Non-ID",
  EFDR_M = NULL,
  EFDR_D = NULL,
  EFDR_R = NULL
)
g1 <- BASiCS_PlotDE(ide1, Plots = c("MAPlot", "VolcanoPlot"))
g2 <- BASiCS_PlotDE(ide2, Plots = c("MAPlot", "VolcanoPlot"))
ggsave(g1, file = "figs/de_id_tung.pdf", width = 12, height = 10)
ggsave(g2, file = "figs/de_id_zeisel.pdf", width = 12, height = 10)



ess_id <- lapply(c("mu", "delta", "epsilon"), 
  function(param) {
    lapply(seq_along(chains_id),
      function(i) {
        x <- chains_id[[i]]
        cbind(
          ess = effectiveSize(x@parameters[[param]]),
          gene = colnames(x@parameters[[param]]),
          param = param,
          desc_id[i, ]
        )
      }
    )
  }
)

ess_dfs <- lapply(ess_id, function(x) do.call(rbind, x))
ess_df_id <- do.call(rbind, ess_dfs)
ess_df_id$data <- gsub("tung", "Tung", ess_df_id$data)
ess_df_id$data <- gsub("zeisel", "Zeisel", ess_df_id$data)
ess_df_id[["param"]] <- factor(ess_df_id[["param"]],
  levels = c("mu", "delta", "epsilon")
)



ggplot(ess_df_id, aes(x = param, y = ess, color = type, fill = type)) + 
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
    breaks = c("id", "non_id"), 
    labels = c("Identifiable", "Non-identifiable"),
    palette = "Dark2"
  )

ggsave("figs/ess_id.pdf", width = 6, height = 4)
