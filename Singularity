Bootstrap: docker
From: rocker/tidyverse:3.6.1

%post
    install2.r --error \
        "devtools" \
        "dplyr" \
        "ggplot2" \
        "here" \
        "magrittr" \
        "RColorBrewer" \
        "readxl" \
        "reshape2" \
        "BiocManager" \
        "viridis" \
        "rstantools" \
        "rstan" \
        "parallelMCMCcombine"
    R -e 'BiocManager::install(
            c(
                "BASiCS",
                "SingleCellExperiment",
                "scater",
                "scran"
            )
        )'
    installGithub.r catavallejos/BASiCS@prior_scaling

%files
    data Scalability/
    data-raw Scalability/
    datasets Scalability/
    DESCRIPTION Scalability/
    .here Scalability/
    inst Scalability/
    man Scalability/
    NAMESPACE Scalability/
    R Scalability/
    src Scalability/
    tests Scalability/
    tools Scalability/



%runscript
    exec /usr/bin/env bash
