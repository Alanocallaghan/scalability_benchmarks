# qsub -t 1-$(cat data/divide_and_conquer_grid.txt | wc -l) -tc 10 scripts/divide_and_conquer.sh
qsub -t 1-$(cat data/advi_grid.txt | wc -l) -tc 5 scripts/advi.sh
