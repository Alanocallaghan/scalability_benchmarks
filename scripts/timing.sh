#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=5G
#$ -cwd
#$ -l h_rt=72:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.5.0

Rscript scripts/timing.R
