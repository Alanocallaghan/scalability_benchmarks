#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=5G
#$ -cwd
#$ -pe sharedmem 5
#$ -l h_rt=72:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/compilers/gcc/5.5.0
module load igmm/apps/R/3.6.1

Rscript scripts/timing.R tung
Rscript scripts/timing.R zeisel
Rscript scripts/timing.R pbmc
Rscript scripts/timing.R buettner
