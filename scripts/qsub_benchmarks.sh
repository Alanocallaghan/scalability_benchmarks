# qsub -t 1-$(cat data/grid_head.tsv | wc -l) scripts/divide_and_conquer.sh
qsub -t 1-$(cat data/datasets.txt | wc -l) scripts/advi.sh
