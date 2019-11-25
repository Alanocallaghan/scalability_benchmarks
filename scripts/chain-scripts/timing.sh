#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -cwd
#$ -l h_rt=100:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.1

# Rscript scripts/chain-scripts/timing.R tung
# Rscript scripts/chain-scripts/timing.R zeisel
Rscript scripts/chain-scripts/timing.R pbmc
# Rscript scripts/chain-scripts/timing.R buettner
