library("SingleCellExperiment")
library("BASiCS")


BASiCS_Draw <- function(
        Chain,
        BatchInfo = gsub(".*_Batch([0-9a-zA-Z])", "\\1", colnames(Chain@parameters[["nu"]])),
        N = sample(nrow(Chain@parameters[["nu"]]), 1)
    ) {
    
    BASiCS_Sim(
        Mu = Chain@parameters[["mu"]][N, ],
        Delta = Chain@parameters[["delta"]][N, ],
        Phi = Chain@parameters[["phi"]][N, ],
        S = Chain@parameters[["s"]][N, ],
        BatchInfo = BatchInfo,
        Theta = Chain@parameters[["theta"]][N, ]
    )
}

density_plot <- function(base_param, ppc_param) {
    df <- do.call(cbind, ppc_param)
    df <- as.data.frame(df)
    df <- cbind(df, base = base_param[[config$data]])
    mdf <- reshape2::melt(as.matrix(df))
    ggplot(mdf, aes(x = value, group = Var2, color = Var2 == "base")) +
        scale_color_brewer(
            palette = "Set1",
            labels = c("Base", "Posterior sample"),
            guide = FALSE
        ) +
        geom_density() +
        scale_x_log10()  
}

## distribution of mean expression
get_means <- function(datas) {
  sapply(datas, function(x) rowMeans(counts(x)), simplify = FALSE)
}

## distribution of proportion of zeros
get_dropout <- function(datas) {
  sapply(datas, function(x) rowMeans(counts(x) != 0), simplify = FALSE)
}

## distribution of cell size
get_size <- function(datas) {
  sapply(datas, function(x) colSums(counts(x)), simplify = FALSE)
}

## distribution of cell complexity
get_complexity <- function(datas) {
  sapply(datas, function(x) colSums(counts(x) != 0), simplify = FALSE)
}



datas <- sapply(
    c("tung", "zeisel", "buettner", "pbmc"),
    function(dataset) {
        readRDS(paste0("data/", dataset, ".rds"))
    },
    simplify = FALSE
)




base_means <- get_means(datas)
base_dropout <- get_dropout(datas)
base_size <- get_size(datas)
base_complexity <- get_complexity(datas)



# files <- list.files("outputs/divide_and_conquer", full.names = TRUE)
# which(sapply(files, function(x) {
#   c <- readRDS(paste0(x, "/config.rds"))
#   c$data == "pbmc" && c$chains == 1
# }))


file <- list.files("outputs/divide_and_conquer", full.names = TRUE)[33]
chains <- readRDS(paste0(file, "/chains.rds"))
config <- readRDS(paste0(file, "/config.rds"))
if (is.list(chains)) {
    chains <- BASiCS:::.combine_subposteriors(
        chains,
        SubsetBy = config$by
    )
}

datas[[config$data]] <- datas[[config$data]][
    intersect(rownames(datas[[config$data]]), rownames(chains)), ]

chains@parameters[c("mu", "delta", "epsilon")] <- lapply(
    chains@parameters[c("mu", "delta", "epsilon")],
    function(x) {
        x[, rownames(datas[[config$data]])]
    }
)


chains@parameters[c("nu", "s")] <- lapply(
    chains@parameters[c("nu", "s")],
    function(x) {
        colnames(x) <- gsub("_Batch.*$", "", colnames(x))
        x
    }
)

datas[[config$data]] <- datas[[config$data]][, 
    intersect(colnames(datas[[config$data]]), colnames(chains@parameters$nu))]

chains@parameters[c("nu", "s")] <- lapply(
    chains@parameters[c("nu", "s")],
    function(x) {
        x[, colnames(datas[[config$data]])]
    }
)


ppc_datas <- lapply(801:850, function(i) {
    suppressMessages(BASiCS_Draw(chains, N = i))
})

ppc_means <- get_means(ppc_datas)
ppc_dropout <- get_dropout(ppc_datas)
ppc_size <- get_size(ppc_datas)
ppc_complexity <- get_complexity(ppc_datas)


g <- density_plot(base_means, ppc_means)
g <- density_plot(base_dropout, ppc_dropout)
g <- density_plot(base_size, ppc_size)
g <- density_plot(base_complexity, ppc_complexity)


# sapply(1:50, 
#     function(i) {
#         ks.test(ppc_means[[i]], base_means[[config$data]])$p.value
#     }
# )

