.PHONY: datasets $(txts)

txts := $(wildcard data/*.txt)

all: outputs

downloads/chen.rds:
	curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/chen.rds > downloads/chen.rds

data/chen.rds: downloads/chen.rds scripts/data-scripts/chen.R
	qsub scripts/data-scripts/data-runner.sh "scripts/data-scripts/chen.R"

data/tung.rds: scripts/data-scripts/tung.R
	qsub scripts/data-scripts/data-runner.sh "scripts/data-scripts/tung.R"

data/zeisel.rds: scripts/data-scripts/zeisel.R
	qsub scripts/data-scripts/data-runner.sh "scripts/data-scripts/zeisel.R"

data/buettner.rds: scripts/data-scripts/buettner.R
	qsub scripts/data-scripts/data-runner.sh "scripts/data-scripts/buettner.R"

# data/pbmc.rds: scripts/data-scripts/pbmc.Rmd
# 	Rscript -e 'rmarkdown::render("scripts/data-scripts/pbmc.Rmd")'

datasets: data/tung.rds data/zeisel.rds data/buettner.rds data/pbmc.rds

data/%.txt: scripts/data-scripts/grids.R
	qsub scripts/data-scripts/data-runner.sh "scripts/data-scripts/grids.R"

outputs/divide_and_conquer: scripts/chain-scripts/divide_and_conquer.sh 
outputs/divide_and_conquer: datasets data/divide_and_conquer_grid.txt
	qsub -t 1-$(shell cat data/divide_and_conquer_grid.txt | wc -l) -tc 10 scripts/chain-scripts/divide_and_conquer.sh

outputs/advi: scripts/chain-scripts/advi.R datasets data/advi_grid.txt
	qsub -t 1-$(shell cat data/advi_grid.txt | wc -l) -tc 5 scripts/chain-scripts/advi.sh 

outputs/downsampling: scripts/chain-scripts/downsampling_divide.sh 
outputs/downsampling: scripts/chain-scripts/downsampling_reference.sh datasets
	qsub -t 1-$(shell cat data/downsampling_grid_ref.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_reference.sh && \
	qsub -t 1-$(shell cat data/downsampling_grid.txt | wc -l) -tc 5 scripts/chain-scripts/downsampling_divide.sh

outputs/identifiable: data/datasets_batch.txt scripts/chain-scripts/identifiable.sh
	qsub -t 1-$(shell cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/identifiable.sh

outputs/batchinfo: data/datasets_batch.txt scripts/chain-scripts/batchinfo.sh
	qsub -t 1-$(shell cat data/datasets_batch.txt | wc -l) -tc 5 scripts/chain-scripts/batchinfo.sh

outputs/time: datasets scripts/chain-scripts/timing.sh
	qsub scripts/chain-scripts/timing.sh

outputs: outputs/divide_and_conquer outputs/advi
outputs: outputs/downsampling outputs/identifiable outputs/batchinfo outputs/time
