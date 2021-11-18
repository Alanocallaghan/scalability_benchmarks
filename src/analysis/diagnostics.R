library("here")
library("BASiCS")

ess <- function (x) {
    vars <- matrixStats::colVars(x)
    spec <- numeric(ncol(x))
    has_var <- vars > 1e-10
    if (any(has_var, na.rm = TRUE)) {
        spec[which(has_var)] <- apply(x[, which(has_var), drop = FALSE],
            2, function(y) {
                a <- ar(y, aic = TRUE)
                a$var.pred / (1 - sum(a$ar))^2
            })
    }
    setNames(ifelse(spec == 0, 0, nrow(x) * vars/spec), colnames(x))
}

plot_all_diags <- function(df, param, measure) {
    measure_name <- BASiCS:::.MeasureName(measure)
    diag_all_list <- parallel::mclapply(
        seq_len(nrow(df)),
        function(i) {
            cat(i, "/", nrow(df), "\n")
            chain <- readRDS(df[[i, "file"]])
            if (length(chain) > 1) {
                suppressMessages(
                    chain <- BASiCS:::.combine_subposteriors(
                        chain,
                        SubsetBy = "gene"
                    )
                )
            }
            data.frame(
                data = df[[i, "data"]],
                chains = df[[i, "chains"]],
                seed = df[[i, "seed"]],
                by = df[[i, "by"]],
                feature = colnames(chain@parameters[[param]]),
                diag = BASiCS:::.GetMeasure(chain, param, measure)
            )
        }, mc.cores = 4
    )
    diag_all <- bind_rows(diag_all_list)
    diag_all[which(diag_all[["chains"]] == 1), "by"] <- "Reference"
    diag_all <- diag_all[diag_all[["by"]] != "advi", ]
    diag_all[diag_all[["by"]] == "gene", "by"] <- "Divide and conquer"
    diag_all[["chains"]] <- factor(
        diag_all[["chains"]],
        levels = sort(unique(diag_all[["chains"]]))
    )

    diag_all[["data"]] <- sub(
        "([\\w])([\\w]+)", "\\U\\1\\L\\2",
        diag_all[["data"]],
        perl = TRUE
    )

    scale <- if (measure == "ess") {
        list(
            scale_y_log10(name = measure_name),
            geom_hline(
                yintercept = 100,
                col = "firebrick",
                linetype = "dashed"
            )
        )
    } else {
        list(
            scale_y_continuous(name = measure_name),
            geom_hline(
                yintercept = c(-3, 3),
                col = "firebrick",
                linetype = "dashed"
            )
        )
    }

    g <- ggplot(diag_all) +
        aes(x = chains, y = diag, color = by, fill = by) +
        # geom_jitter(height = 0, width = 0.2) +
        geom_violin(alpha = 0.2) +
        geom_boxplot(alpha = 0.2, width = 0.1, outlier.colour = NA) +
        facet_wrap(~data, scales = "free_y") +
        scale_fill_brewer(
            name = "Inference method",
            palette = "Dark2",
            aesthetics = c("fill", "color")
        ) +
        scale +
        theme(panel.grid = element_blank()) +
        labs(
            x = "Partitions",
            y = "Effective sample size"
        )
    dir.create(
        sprintf("figs/%s", gsub("\\.", "_", measure)),
        showWarnings = FALSE
    )
    ggsave(g,
        file = sprintf("figs/%s/%s_all.pdf", gsub("\\.", "_", measure), param),
        width = 5, height = 4
    )
    invisible(g)
}

plot_all_diags(df, "mu", "ess")
plot_all_diags(df, "delta", "ess")
plot_all_diags(df, "epsilon", "ess")
plot_all_diags(df, "mu", "geweke.diag")
plot_all_diags(df, "delta", "geweke.diag")
plot_all_diags(df, "epsilon", "geweke.diag")
