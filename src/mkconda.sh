# grep -rh "library(" | sed -e 's/^[ \t]*//' | sort | uniq
# library("argparse")
# library("coda")
# library("dplyr")
# library("ggbeeswarm")
# library("ggplot2")
# library("here")
# library("viridis")
# library("BASiCStan")
# library("BASiCS")
# library("BiocParallel")
# library("scater")
# library("scran")
# library("SingleCellExperiment")

if [ -f "/etc/profile.d/modules.sh" ]; then
    module load roslin/gcc/7.3.0
fi

set -eu                                                                                                                                                                   
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then                                                                                                                       
    source ~/miniconda3/etc/profile.d/conda.sh                                                                                                                            
fi                                                                                                                                                                        
set +eu

conda create -y -n scalability \
    r-base=4.3.1 \
    r-argparse \
    r-coda \
    r-dplyr \
    r-ggbeeswarm \
    r-ggplot2 \
    r-here \
    r-viridis \
    snakemake \
    r-biocmanager \
    r-devtools \
    r-rstan \
    r-ggpointdensity \
    r-rcpparmadillo \
    r-ggrastr \
    # bioconductor-scrnaseq
    # bioconductor-basics \
    # bioconductor-scater \
    # bioconductor-scran \
    # bioconductor-singlecellexperiment \

conda init
conda activate scalability
Rscript -e 'BiocManager::install(c("BASiCS", "BASiCStan", "BiocParallel", "scater", "scran", "SingleCellExperiment", "scRNAseq"), version=3.18)'
