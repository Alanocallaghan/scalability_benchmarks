library("here")
source(here("src/analysis/preamble.R"))

time_files <- list.files(here("outputs/time/"), full.names = TRUE)

time_data <- do.call(c,
    sapply(time_files,
    function(x) {
        # print(x)
        mean(readRDS(x))
    },
    simplify = FALSE)
)
time_df_dc <- data.frame(
    data = gsub(
        ".*(chen|buettner|zeisel|tung|ibarra-soria).*",
        "\\1",
        names(time_data)
    ),
    chains = gsub(
        ".*(chen|buettner|zeisel|tung|ibarra-soria)_(\\d+).rds", "\\2",
        names(time_data)
    ),
    time = time_data,
    row.names = NULL
)
time_df_dc <- time_df_dc[time_df_dc$data != "ibarra-soria", ]
time_df_dc <- time_df_dc[time_df_dc$chains != 128, ]

advi_time_files <- list.files(
    "outputs/advi",
    recursive = TRUE,
    pattern = "time.rds",
    full.names = TRUE
)
                                                                 
advi_time_df <- data.frame(
    data = gsub(
        ".*/data-(\\w+)_seed-(\\d+)/time.rds", "\\1",
        advi_time_files
    ),
    seed = gsub(
        ".*/data-(\\w+)_seed-(\\d+)/time.rds", "\\2",
        advi_time_files
    ),
    file = advi_time_files
)
advi_time_df$time <- sapply(advi_time_df$file,
    function(x) {
        readRDS(x)[["elapsed"]]
    }
)

advi_time_df <- merge(advi_time_df, data_dims)

advi_time_df <- advi_time_df %>%
    dplyr::group_by(data) %>%
    dplyr::summarise(
        .groups = "drop_last",
        time = median(time),
        nGenes = nGenes[[1]],
        nCells = nCells[[1]],
    )

advi_time_df$data <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    advi_time_df$data,
    perl = TRUE
)


time_df_merge <- merge(time_df_dc, data_dims, all = TRUE)
time_df_merge <- time_df_merge[, 
    c("data", "chains", "time", "nGenes", "nCells")
]
time_df_merge <- time_df_merge %>%
    group_by(data) %>%
    mutate(
        nGenes = mean(nGenes, na.rm = TRUE),
        nCells = mean(nCells, na.rm = TRUE)
    )

time_df_merge$data <- sub(
    "([\\w])([\\w]+)", "\\U\\1\\L\\2",
    time_df_merge$data,
    perl = TRUE
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
            color = data
            # color = paste(
            #   data, "\n", 
            #   nGenes, "genes;", nCells, "cells",
            #   "\n"
            # )
        )
    ) +
    geom_line(aes(group = data)) +
    geom_hline(
        aes(
        yintercept = time / 60,
        color =data
        # color = paste(
        #   data, "\n", 
        #   nGenes, "genes;",
        #   nCells, "cells",
        #   "\n"
        # ),
        ),
        linetype = "dashed",
        data = advi_time_df
    ) +
    scale_x_continuous(
        name = "Number of partitions",
        trans = "log2",
        breaks = c(1, 2, 4, 8, 16, 32, 64, 128)
    ) +
    scale_y_continuous(name = "Time (mins)", trans = "log2") +
    scale_color_brewer(name = "Data", palette = "Set2") +
    # scale_linetype_discrete(
    #   name = "Method", limits = c("Divide and\nconquer", "ADVI")
    # ) +
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom"
    )
    #  +
    # guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    # theme(
    #   legend.position = "bottom",
    #   legend.box = "vertical"
    # )

ggsave(
    file = here("figs/time_plot.pdf"),
    width = 4,
    height = 4
)
