qsub -t 1-$(cat data/divide_and_conquer_grid.tx | wc -l) scripts/divide_and_conquer.sh
qsub -t 1-$(cat data/advi_grid.txt | wc -l) scripts/advi.sh
