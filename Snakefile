chains = [1, 2, 4, 8, 16, 32, 64]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
# seeds = [14, 21, 28, 35, 42]
seeds = [42]
fractions = [x/10 for x in range(2, 11, 2)]
iterations = 40

# configfile: "config/snakemake_config.yaml"
# conda: config["conda"]

shell.prefix("source src/modules.sh; ")

rule all:
    input:
        expand(
            "outputs/divide_and_conquer/data-{dataset}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/",
            dataset = data,
            nsubsets = chains,
            by = by,
            seed = seeds
        ),
        # expand(
        #     "outputs/advi/data-{dataset}_seed-{seed}/",
        #     data = data,
        #     seed  = seeds
        # ),
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
        )


# rule plots: ## todo

rule advi:
    conda:
         "conda.yaml"
    output:
        "outputs/advi/data-{dataset}_seed-{seed}/"
    input:
        "rdata/{dataset}.rds"
    shell:
        """
        Rscript ./src/chain/advi.R \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --output {output}
        """

rule time:
    # conda: "conda.yaml"
    output:
        expand(
            "outputs/time/{dataset}_{n}.rds",
            dataset = data
            n = [2, 4, 8, 16, 32]
        )
    input:
        "rdata/{dataset}.rds"
    shell:
        """
        Rscript ./src/chain
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
        Rscript ./src/chain/divide_and_conquer.R \
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
        Rscript ./src/chain/downsampling_reference.R \
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
        Rscript ./src/chain/downsampling_divide.R \
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
        Rscript ./src/chain/removing_cells_ref.R \
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
        Rscript ./src/chain/removing_cells.R \
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
        Rscript ./src/chain/true_positives.R \
            --iterations {iterations} \
            --nsubsets {wildcards.nsubsets} \
            --seed {wildcards.seed} \
            --output {output}
        """


rule batchinfo:
    # conda:
    #     "conda.yaml"
    resources: mem_mb=10000, runtime=3000
    input:
        "rdata/{dataset}.rds"
    output:
        directory("outputs/batchinfo/data-{dataset}/")
    shell:
        """
        Rscript ./src/chain/batchinfo.R \
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
        Rscript ./src/chain/batchinfo.R \
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

