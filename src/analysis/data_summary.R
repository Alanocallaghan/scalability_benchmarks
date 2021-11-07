library("SingleCellExperiment")
library("argparse")
library("xtable")

parser <- ArgumentParser()
parser$add_argument("-i", "--input")

args <- parser$parse_args()

sces <- args[["input"]]
if (is.null(sces)) {
    sces <- list.files("rdata", full.names = TRUE)
}

dims <- lapply(sces, function(x) dim(readRDS(x)))

df <- data.frame(Dataset = gsub("rdata/", "", sces))
df$Dataset <- gsub(".rds", "", df$Dataset)
df$Dataset <- gsub("([\\w])([\\w]+)", "\\U\\1\\L\\2", df$Dataset, perl = TRUE)
df <- cbind(df, as.data.frame(do.call(rbind, dims)))
colnames(df)[2:3] <- c("Genes", "Cells")

out <- xtable(
    df,
    align = "rr|rr",
    caption = "The number of cells and genes present in each dataset following the
    filtering steps described in \\cref{sec:DataProcessing}.",
    label = "tab:CellsGenes"
)
print(out,
    include.rownames = FALSE,
    file = "tables/data-summary.tex"
)
