#!/bin/bash -f
#$ -S /bin/bash
#$ -pe sharedmem 16
#$ -l h_vmem=5G
#$ -cwd
#$ -l h_rt=24:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.1

Rscript scripts/chain-scripts/cell_splitting.R
