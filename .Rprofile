if ("BiocParallel" %in% rownames(utils::installed.packages())) {
    library("BiocParallel")
    register(SerialParam())
}
options(menu.graphics=FALSE)
options(mc.cores=1)
options(Ncpus=4)
## Github personal access token
source("PAT.R")
