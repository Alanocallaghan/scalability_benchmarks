# grep -rh "library(" | sed -e 's/^[ \t]*//' | sort | uniq
# library("argparse")
# library("BASiCS")
# library("coda")
# library("dplyr")
# library("edgeR")
# library("ggbeeswarm")
# library("ggplot2")
# library("here")
# library("Scalability")
# library("scater")
# library("scran")
# library("Seurat")
# library("SingleCellExperiment")
# library("viridis")

mamba create -n scalability r-argparse \
    bioconductor-basics \
    r-coda \
    r-dplyr \
    bioconductor-edgeR \
    r-ggbeeswarm \
    r-ggplot2 \
    r-here \
    bioconductor-scater \
    bioconductor-scran \
    bioconductor-singlecellexperiment \
    r-viridis \
    snakemake \
    r-biocmanager \
    r-devtools \
    r-rstan

## also need to devtools::install_github("Alanocallaghan/BASiCStan")
conda activate scalability
Rscript -e 'devtools::install_github("Alanocallaghan/BASiCStan")'
Rscript -e 'devtools::install_github("catavallejos/BASiCS", ref="hotfix")'

