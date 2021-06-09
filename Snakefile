chains = [1, 2, 4, 8, 16, 32, 64]
chains_timing = [1, 2, 4, 8, 16, 32, 64, 128]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_spikes = ["buettner", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
seeds = [14, 21, 28, 35, 42]
# seeds = [42]
fractions = [x/10 for x in range(2, 11, 2)]
iterations = 40000

# configfile: "config/snakemake_config.yaml"
# conda: config["conda"]

shell.prefix("source src/modules.sh; ")

rule all:
    input:
        "figs/diffexp_plot.pdf",
        "figs/overlap_diff_genes.pdf",
        "figs/libsize_density.pdf",
        "figs/complexity_density.pdf",
        "figs/expression_density.pdf",
        "figs/dropout_density.pdf",
        "figs/true-positives.pdf",
        "figs/time_plot.pdf",
        "figs/de_batch_tung.pdf",
        "figs/de_batch_zeisel.pdf",
        "figs/ess_batch.pdf",
        "figs/removing_cells.pdf",
        "figs/de_id_tung.pdf",
        "figs/de_id_zeisel.pdf",
        "figs/ess_id.pdf",
        "figs/downsampling.pdf"


rule plots: ## todo
    resources: mem_mb=50000, runtime=10000
    input:
        expand(
            "outputs/divide_and_conquer/data-{dataset}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/",
            dataset = data,
            nsubsets = chains,
            by = by,
            seed = seeds
        ),
        expand(
            "outputs/advi/data-{dataset}_seed-{seed}/",
            dataset = data_spikes,
            seed  = seeds
        ),
        expand(
            "outputs/downsampling/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/",
            dataset = data,
            seed = seeds,
            fraction = fractions
        ),
        expand(
            "outputs/downsampling/reference/data-{dataset}_fraction-{fraction}/",
            dataset = data,
            fraction = fractions
        ),
        expand(
            "outputs/removing/reference/data-{dataset}_fraction-{fraction}/",
            dataset = data,
            fraction = fractions
        ),
        expand(
            "outputs/removing/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/",
            dataset = data,
            seed = seeds,
            fraction = fractions
        ),
        expand(
            "outputs/batchinfo/data-{dataset}/",
            dataset = data_batch
        ),
        expand(
            "outputs/true-positives/data-ibarra-soria_nsubsets-{nsubsets}_seed-{seed}.rds",
            nsubsets = chains,
            seed = seeds
        ),
        expand(
            "outputs/time/{dataset}_{n}.rds",
            dataset = data,
            n = chains_timing
        )
    output:
        "figs/diffexp_plot.pdf",
        "figs/overlap_diff_genes.pdf",
        "figs/libsize_density.pdf",
        "figs/complexity_density.pdf",
        "figs/expression_density.pdf",
        "figs/dropout_density.pdf",
        "figs/true-positives.pdf",
        "figs/time_plot.pdf",
        "figs/de_batch_tung.pdf",
        "figs/de_batch_zeisel.pdf",
        "figs/ess_batch.pdf",
        "figs/removing_cells.pdf",
        "figs/de_id_tung.pdf",
        "figs/de_id_zeisel.pdf",
        "figs/ess_id.pdf",
        "figs/downsampling.pdf"
    shell:
        """
        Rscript src/analysis/main.R
        """

rule advi:
    conda:
         "conda.yaml"
    resources: mem_mb=20000
    output:
        directory("outputs/advi/data-{dataset}_seed-{seed}/")
    input:
        "rdata/{dataset}.rds"
    shell:
        """
        Rscript ./src/chains/advi.R \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --output {output}
        """

rule time:
    # conda: "conda.yaml"
    resources: mem_mb=10000, runtime=10000
    output:
        "outputs/time/{dataset}_{n}.rds"
    input:
        "rdata/{dataset}.rds"
    shell:
        """
        Rscript ./src/chains/timing.R \
            --data {wildcards.dataset} \
            --nsubsets {wildcards.n} \
            --iterations {iterations} \
            --output {output}
        """


rule divide_and_conquer:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/divide_and_conquer/data-{dataset}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/")
    shell:
        """
        Rscript ./src/chains/divide_and_conquer.R \
            --iterations {iterations} \
            --data {wildcards.dataset} \
            --nsubsets {wildcards.nsubsets} \
            --seed {wildcards.seed} \
            --subsetby {wildcards.by} \
            --output {output}
        """

rule downsampling_ref:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/downsampling/reference/data-{dataset}_fraction-{fraction}/")
    shell:
        """
        Rscript ./src/chains/downsampling_reference.R \
            --iterations {iterations} \
            --data {wildcards.dataset} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule downsampling_divide:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/downsampling/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/")
    shell:
        """
        Rscript ./src/chains/downsampling_divide.R \
            --iterations {iterations} \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule removing_ref:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/removing/reference/data-{dataset}_fraction-{fraction}/")
    shell:
        """
        Rscript ./src/chains/removing_cells_ref.R \
            --iterations {iterations} \
            --data {wildcards.dataset} \
            --fraction {wildcards.fraction} \
            --output {output}
        """

rule removing_divide:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/removing/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/")
    shell:
        """
        Rscript ./src/chains/removing_cells.R \
            --iterations {iterations} \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """



rule true_positives:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=5000
    input:
        "rdata/ibarra-soria.rds"
    output:
        "outputs/true-positives/data-ibarra-soria_nsubsets-{nsubsets}_seed-{seed}.rds"
    shell:
        """
        Rscript ./src/chains/true_positives.R \
            --iterations {iterations} \
            --nsubsets {wildcards.nsubsets} \
            --seed {wildcards.seed} \
            --output {output}
        """


rule batchinfo:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=20000, runtime=5000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/batchinfo/data-{dataset}/")
    shell:
        """
        Rscript ./src/chains/batchinfo.R \
	    --iterations {iterations} \
            --data {wildcards.dataset} \
            --output {output}
        """


rule cell:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=3000
    input:
        "rdata/{dataset}.rds"
    output:
        "outputs/cell_splitting/{dataset}.rds"
    shell:
        """
        Rscript ./src/chains/batchinfo.R \
	    --iterations {iterations} \
            --data {wildcards.dataset} \
            --output {output}
        """


rule data:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "src/data/{dataset}.R"
    output:
        "rdata/{dataset}.rds"
    shell:
        """
        Rscript {input}
        """

