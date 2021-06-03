time_files <- list.files(here("outputs/time/"), full.names = TRUE)

time_df_dc <- data.frame(
  data = gsub(".*/(\\w+)_(\\d+).rds", "\\1", time_files),
  chains = gsub(".*/(\\w+)_(\\d+).rds", "\\2", time_files)
)
time_df_dc <- time_df_dc[rep(seq_len(nrow(time_df_dc)), each = 6), ]
time_df_dc[["times"]] <- as.vector(sapply(time_files, readRDS))
time_df_dc[["seeds"]] <- seq(7, 42, length.out = 6)



advi_time_df <- df %>% 
  dplyr::filter(by == "advi") %>%
  dplyr::group_by(data) %>% 
  dplyr::summarise(
    .groups = "drop_last",
    time = median(time),
    nGenes = nGenes[[1]],
    nCells = nCells[[1]],
  )

advi_time_df$data <- sub(
  "([[:alpha:]])", "\\U\\1",
  advi_time_df$data,
  perl = TRUE
)
advi_time_df$data <- sub(
  "Pbmc",
  "10x PBMC",
  advi_time_df$data
)


time_df <- df %>% 
  dplyr::filter(by != "advi") %>%
  dplyr::group_by(data, chains) %>% 
  dplyr::summarise(
    .groups = "drop_last",
    time = median(time),
    nGenes = nGenes[[1]],
    nCells = nCells[[1]],
  )

time_df_merge <- merge(time_df, time_df_dc, all = TRUE)
time_df_merge <- time_df_merge[, 
  c("data", "chains", "seeds", "time", "times", "nGenes", "nCells")
]
time_df_merge <- time_df_merge %>% 
  group_by(data) %>%
  mutate(
    nGenes = mean(nGenes, na.rm = TRUE),
    nCells = mean(nCells, na.rm = TRUE)
  )
ind_na_time <- is.na(time_df_merge[["times"]])
time_df_merge[ind_na_time, "times"] <- time_df_merge[ind_na_time, "time"]
time_df_merge[["time"]] <- time_df_merge[["times"]]

time_df_merge$data <- sub(
  "([[:alpha:]])", "\\U\\1",
  time_df_merge$data,
  perl = TRUE
)
time_df_merge$data <- sub(
  "Pbmc",
  "10x PBMC",
  time_df_merge$data
)

time_df <- time_df_merge %>% 
  dplyr::group_by(data, chains) %>% 
  dplyr::summarise(
    .groups = "drop_last",
    time = median(time),
    nGenes = nGenes[[1]],
    nCells = nCells[[1]],
  )


g <- ggplot(
  time_df, 
  aes(
    x = as.numeric(chains),
    y = time / 60,
    color = paste(
      data, "\n", 
      nGenes, "genes;", nCells, "cells",
      "\n"
    )
  )
) +
  geom_line(aes(group = data, linetype = "Divide and\nconquer")) +
  geom_hline(
    aes(
      yintercept = time / 60,
      color = paste(
        data, "\n", 
        nGenes, "genes;",
        nCells, "cells",
        "\n"
      ),
      lty = "ADVI",
    ),
    data = advi_time_df
  ) +
  scale_x_continuous(
    name = "Partitions",
    trans = "log2",
    breaks = c(1, 2, 4, 8, 16, 32)
  ) +
  scale_y_continuous(name = "Time (mins)", trans = "log10") +
  scale_color_brewer(name = "Data", palette = "Set2") +
  scale_linetype_discrete(name = "Method", limits = c("Divide and\nconquer", "ADVI"))

ggsave(
  file = here("figs/time_plot.pdf"),
  width = 7,
  height = 7
)



# mean_df <- df %>% dplyr::filter(by != "advi") %>%
#   dplyr::group_by(data, chains) %>% 
#   dplyr::summarize(
#     time = mean(time),
#     nGenes = nGenes[[1]],
#     nCells = nCells[[1]],
#   )


# df$chains <- factor(df$chains, levels = c(1, 2, 4, 8, 16, 32))

# g <- ggplot(
#   df,
#   aes(
#     x = chains,
#     y = time / 3600,
#     color = paste(
#       data, "\n", 
#       nGenes, "genes;", nCells, "cells",
#       "\n"
#     )
#   )
# ) +
#   geom_line(data = mean_df, aes(group = data, linetype = "Divide and\nconquer")) +
#   geom_boxplot()
#   geom_hline(
#     aes(
#       yintercept = time / 3600,
#       color = paste(
#         data, "\n", 
#         nGenes, "genes;",
#         nCells, "cells",
#         "\n"
#       ),
#       lty = "ADVI",
#     ),
#     data = advi_time_df
#   ) +
#   scale_x_continuous(
#     name = "Partitions",
#     trans = "log2",
#     breaks = c(1, 2, 4, 8, 16, 32)
#   ) +
#   scale_y_continuous(name = "Time (hr)") +
#   scale_color_brewer(name = "Data", palette = "Set2") +
#   scale_linetype_discrete(name = "Method", limits = c("Divide and\nconquer", "ADVI"))

# g <- ggplot(out_f, 
#     aes(
#       x = chains, 
#       y = time / 3600, 
#     )
#   ) +
#   geom_point() +
#   # geom_hline(yintercept = 12255 / 3600, lty = "dashed", colour = "grey60") + 
#   # annotate(x = 1, y = 12255 / 3600 * 1.2, size = fs,
#   #   label = "Time taken for ADVI", geom = "text", hjust = 0) +
#   geom_hline(yintercept = tmax / (10 * 3600), lty = "dashed", colour = "grey60") + 
#   annotate(x = 1, y = tmax / (10 * 3600) * 1.2, size = fs,
#     label = "10x speedup", geom = "text", hjust = 0) +
#   geom_hline(yintercept = tmax / (100 * 3600), lty = "dashed", colour = "grey60") + 
#   annotate(x = 1, y = tmax / (100 * 3600) * 1.2, size = fs,
#     label = "100x speedup", geom = "text", hjust = 0) +
#   geom_line() +
#   geom_text(
#     data = out_f, 
#     hjust = 0.4,
#     size = fs,
#     mapping = aes(
#       x = chains,
#       y = 25,
#       label = n_str),
#     show.legend = FALSE
#   ) +
#   labs(x = "Number of partitions", y = "Time") +
#   scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64, 128)) +
#   scale_y_continuous(
#     trans = "log10", 
#     breaks = c(0.16666, 0.5, 1, 2, 5, 15), 
#     labels = c("10 min", "30 min", "1 hr", "2 hr", "5 hr", "15 hr")
#   )
