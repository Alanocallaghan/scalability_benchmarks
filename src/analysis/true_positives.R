library("BASiCS")
library("ggplot2")
library("here")
source(here("src/analysis/functions.R"))

files <- list.files(
    "outputs/true-positives/",
    pattern = ".*.rds",
    full.names = TRUE
)

data <- lapply(files, readRDS)

pos_metadata <- data.frame(
    data = gsub(".*data-(\\w+-\\w+).*", "\\1", files),
    nsubsets = gsub(".*nsubsets-(\\d+).*", "\\1", files),
    seed = gsub(".*seed-(\\d+).*", "\\1", files)
)

pos_metadata$test <- lapply(data, function(x) x$test)
params <- c("Mean", "Disp", "ResDisp")
pg <- paste0(params, "Gene")
pos_metadata[, pg] <- NA
for (i in 1:nrow(pos_metadata)) {
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
        pos_metadata[i, pg[[j]]] <- as.data.frame(df[[j]])
    }
}

ref <- pos_metadata[pos_metadata$nsubsets == 1, ]
jp <- paste0(params, "Jaccard")
pos_metadata[, jp] <- NA
for (i in 1:nrow(pos_metadata)) {
    for (j in seq_along(jp)) {
        pos_metadata[i, jp[[j]]] <- jaccard(
            ref[1, pg[[j]]],
            pos_metadata[i, pg[[j]]]
        )
    }
}

jd <- pos_metadata[, c("nsubsets", jp)]
mdf <- reshape2::melt(jd, id.var="nsubsets")
mdf$variable <- gsub("Jaccard", "", mdf$variable)
mdf$nsubsets <- factor(mdf$nsubsets, levels = sort(unique(as.numeric(mdf$nsubsets))))

g <- ggplot(mdf) +
    aes(nsubsets, y = value, colour = variable) +
    geom_point() +
    labs(x = "Number of subsets", y = "Jaccard Index") +
    ylim(0, 1) +
    scale_color_brewer(palette="Set1", name = "Parameter") +
    theme_bw()

ggsave(file="figs/true-positives.pdf", width = 7, height = 7)
