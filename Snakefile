chains = [1, 2, 4, 8, 16, 32, 64, 128]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
seeds = [7, 14, 21, 28, 35, 42]
fractions = [x/10 for x in range(2, 11, 2)]

# configfile: "config/snakemake_config.yaml"
# conda: config["conda"]
conda: "env.yaml"

rule all:
    input:
        expand(
            "outputs/divide_and_conquer/data-{data}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/",
            data = data,
            nsubsets = chains,
            by = by,
            seed = seeds
        ),
        # expand(
        #     "outputs/advi/data-{data}_seed-{seed}/",
        #     data = data,
        #     seed  = seeds
        # ),
        expand(
            "outputs/downsampling/divide/data-{data}_fraction-{fraction}_seed-{seed}/",
            data = data,
            seed = seeds,
            fraction = fractions
        ),
        expand(
            "outputs/downsampling/reference/data-{data}_fraction-{fraction}/",
            data = data,
            fraction = fractions
        ),
        expand(
            "outputs/removing/reference/data-{data}_fraction-{fraction}/",
            data = data,
            fraction = fractions
        ),
        expand(
            "outputs/removing/divide/data-{data}_fraction-{fraction}_seed-{seed}/",
            data = data,
            seed = seeds,
            fraction = fractions
        ),
        expand(
            "outputs/batchinfo/data-{data}/",
            data = data_batch
        )

# rule advi:
#     output:
#         "outputs/advi/data-{data}_seed-{seed}/"
#     shell:
#         """
#         ./src/chain-scripts/advi.R \
#             --data {wildcards.data} \
#             --seed {wildcards.seed} \
#             --output {output}
#         """


rule divide_and_conquer:
    # threads: 4
    # resources:
    #     mem_mb=5000
    output:
        "outputs/divide_and_conquer/data-{data}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/"
    shell:
        """
        ./src/chain-scripts/divide_and_conquer.R \
            --data {wildcards.data} \
            --nsubsets {wildcards.nsubsets} \
            --seed {wildcards.seed} \
            --subsetby {wildcards.by} \
            --output {output}
        """

rule downsampling_ref:
    output:
        "outputs/downsampling/reference/data-{data}_fraction-{fraction}/"
    shell:
        """
        ./src/chain-scripts/downsampling_reference.R \
            --data {wildcards.data} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule downsampling_divide:
    output:
        "outputs/downsampling/divide/data-{data}_fraction-{fraction}_seed-{seed}/"
    shell:
        """
        ./src/chain-scripts/downsampling_divide.R \
            --data {wildcards.data} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """


rule removing_ref:
    output:
        "outputs/removing/reference/data-{data}_fraction-{fraction}/"
    shell:
        """
        ./src/chain-scripts/removing_cells_ref.R \
            --data {wildcards.data} \
            --fraction {wildcards.fraction} \
            --output {output}
        """

rule removing_divide:
    output:
        "outputs/removing/divide/data-{data}_fraction-{fraction}_seed-{seed}/"
    shell:
        """
        ./src/chain-scripts/removing_cells.R \
            --data {wildcards.data} \
            --seed {wildcards.seed} \
            --fraction {wildcards.fraction} \
            --output {output}
        """

rule batchinfo:
    output:
        "outputs/batchinfo/data-{data}/"
    shell:
        """
        ./src/chain-scripts/batchinfo.R \
            --data {wildcards.data} \
            --output {output}
        """

rule data:
    input:
        "scripts/data-scripts/{dataset}.R"
    output:
        "data/{dataset}.rds"
    shell:
        """
        Rscript {input}
        """

