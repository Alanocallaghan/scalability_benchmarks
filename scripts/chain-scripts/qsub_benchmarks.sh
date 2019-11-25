# qsub -t 1-$(cat data/divide_and_conquer_grid.txt | wc -l) -tc 10 scripts/chain-scripts/divide_and_conquer.sh
# qsub -t 1-$(cat data/advi_grid.txt | wc -l) -tc 5 scripts/chain-scripts/advi.sh
qsub -t 1-$(cat data/downsampling_grid_ref.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_reference.sh
qsub -t 1-$(cat data/downsampling_grid.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_divide.sh
qsub -t 1-$(cat data/removing_grid_ref.txt | wc -l) -tc 5 scripts/chain-scripts/removing_cells_ref.sh
qsub -t 1-$(cat data/removing_grid.txt | wc -l) -tc 5 scripts/chain-scripts/removing_cells.sh
# qsub -t 1-$(cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/identifiable.sh
# qsub -t 1-$(cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/batchinfo.sh
# qsub scripts/chain-scripts/timing.sh
