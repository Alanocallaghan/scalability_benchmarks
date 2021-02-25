#!/bin/sh

if [ -f "/etc/profile.d/modules.sh" ]; then

    . /etc/profile.d/modules.sh

    module load phys/compilers/gcc/10.2.0

    export C_INCLUDE_PATH=/exports/igmm/software/pkg/el7/apps/hdf5/1.8.13/include

    set +eu

    if [ -f "/home/s1372510/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/s1372510/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/s1372510/miniconda3/bin:$PATH"
    fi

    set -eu

fi


conda activate $(find .snakemake/conda/ -mindepth 1 -maxdepth 1 -type d)

