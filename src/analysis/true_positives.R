library("BASiCS")
library("ggplot2")
library("here")
source(here("src/analysis/functions.R"))

files <- list.files(
    "outputs/true-positives/",
    pattern = ".*.rds",
    full.names = TRUE
)

pos_metadata <- data.frame(
    data = gsub(".*data-(\\w+-\\w+).*", "\\1", files),
    nsubsets = gsub(".*nsubsets-(\\d+).*", "\\1", files),
    seed = gsub(".*seed-(\\d+).*", "\\1", files)
)

pos_metadata$test <- lapply(seq_along(files),
    function(i) {
        cat(i, "/", length(files), "\n")
        readRDS(files[[i]])$test
    }
)
params <- c("Mean", "Disp", "ResDisp")
pg <- paste0(params, "Gene")
pos_metadata[, pg] <- NA
pos_metadata[, pg] <- lapply(pos_metadata[, pg], function(x) rep(list(), length(x)))

for (i in 1:nrow(pos_metadata)) {
    cat(i, "/", nrow(pos_metadata), "\n")
    df <- lapply(
        params,
        function(x) {
            out <- as.data.frame(pos_metadata$test[[i]], Parameter=x)$GeneName
            if (length(out)) {
                out
            } else {
                NA
            }
        }
    )
    for (j in seq_along(pg)) {
        pos_metadata[[i, pg[[j]]]] <- list(df[[j]])
    }
}

ref <- pos_metadata[pos_metadata$nsubsets == 1, ]
pos_metadata_test <- pos_metadata[pos_metadata$nsubsets != 1, ]
jp <- paste0(params, "Jaccard")
pos_metadata_test[, jp] <- NA
for (i in 1:nrow(pos_metadata_test)) {
    cat(i, "/", nrow(pos_metadata_test), "\n")
    for (j in seq_along(jp)) {
        pos_metadata_test[i, jp[[j]]] <- jaccard(
            ref[[1, pg[[j]]]][[1]],
            pos_metadata_test[[i, pg[[j]]]][[1]]
        )
    }
}

jd <- pos_metadata_test[, c("nsubsets", jp)]
mdf <- reshape2::melt(jd, id.var="nsubsets")
mdf$variable <- gsub("Jaccard", "", mdf$variable)
mdf$nsubsets <- factor(mdf$nsubsets, levels = sort(unique(as.numeric(mdf$nsubsets))))

mdf$variable <- plyr::revalue(
    mdf$variable,
    replace = c(
        "Mean" = "mu",
        "Disp" = "delta",
        "ResDisp" = "epsilon"
    )
)
g <- ggplot(mdf) +
    aes(nsubsets, y = value, colour = variable) +
    geom_quasirandom(dodge.width = 0.3, size = 0.8, width = 0.05, groupOnX = TRUE) +
    labs(x = "Number of subsets", y = "Jaccard Index") +
    ylim(0, 1) +
    scale_color_brewer(palette="Set1", name = "Parameter") +
    theme_bw()

ggsave(file="figs/true_positives.pdf", width = 7, height = 7)
