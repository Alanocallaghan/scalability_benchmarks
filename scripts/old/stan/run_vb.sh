#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -pe sharedmem 6
#$ -cwd
#$ -l h_rt=200:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.5.0

Rscript ./data-raw/stan/vb_testde.R
