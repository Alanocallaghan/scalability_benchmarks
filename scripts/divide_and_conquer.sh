#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=6G
#$ -pe sharedmem 16
#$ -cwd
#$ -l h_rt=72:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.0

settings=($(sed -n "$SGE_TASK_ID p" data/divide_and_conquer_grid.txt))

Rscript scripts/divide_and_conquer.R ${settings[@]} /exports/eddie/scratch/s1372510/divide_and_conquer/$SGE_TASK_ID


