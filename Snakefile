chains = [1, 2, 4, 8, 16, 32, 64, 128]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
seeds = [7, 14, 21, 28, 35, 42]
fractions = [x/10 for x in range(2, 11, 2)]

# configfile: "config/snakemake_config.yaml"
# conda: config["conda"]

shell.prefix("source src/modules.sh;")

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
        )

# rule advi:
#     conda:
#          "conda.yaml"
#     output:
#         "outputs/advi/data-{dataset}_seed-{seed}/"
#     input:
#         "data/{dataset}.rds"
#     shell:
#         """
#         Rscript ./src/chain-scripts/advi.R \
#             --data {wildcards.dataset} \
#             --seed {wildcards.seed} \
#             --output {output}
#         """


rule divide_and_conquer:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/divide_and_conquer/data-{dataset}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/")
    shell:
        """
        Rscript ./src/chain-scripts/divide_and_conquer.R \
            --data {wildcards.dataset} \
            --nsubsets {wildcards.nsubsets} \
            --seed {wildcards.seed} \
            --subsetby {wildcards.by} \
            --output {output}
        """

rule downsampling_ref:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/downsampling/reference/data-{dataset}_fraction-{fraction}/")
    shell:
        """
        Rscript ./src/chain-scripts/downsampling_reference.R \
            --data {wildcards.dataset} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule downsampling_divide:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/downsampling/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/")
    shell:
        """
        Rscript ./src/chain-scripts/downsampling_divide.R \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule removing_ref:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/removing/reference/data-{dataset}_fraction-{fraction}/")
    shell:
        """
        Rscript ./src/chain-scripts/removing_cells_ref.R \
            --data {wildcards.dataset} \
            --fraction {wildcards.fraction} \
            --output {output}
        """

rule removing_divide:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/removing/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/")
    shell:
        """
        Rscript ./src/chain-scripts/removing_cells.R \
            --data {wildcards.dataset} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """

rule batchinfo:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "data/{dataset}.rds"
    output:
        directory("outputs/batchinfo/data-{dataset}/")
    shell:
        """
        Rscript ./src/chain-scripts/batchinfo.R \
            --data {wildcards.dataset} \
            --output {output}
        """

rule data:
    # conda:
    #     "conda.yaml"
    resources:
        mem_mb=10000
    input:
        "src/data-scripts/{dataset}.R"
    output:
        "data/{dataset}.rds"
    shell:
        """
        Rscript {input}
        """

