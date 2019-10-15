#!/bin/bash -f
#$ -S /bin/bash
#$ -l h_vmem=30G
#$ -cwd
#$ -l h_rt=72:00:00

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load igmm/apps/singularity/2.6.0

dataset=$(sed -n "$SGE_TASK_ID p" data/datasets.txt)

singularity exec image.sif Rscript scripts/advi.R $dataset /exports/eddie/scratch/s1372510/advi/$SGE_TASK_ID
