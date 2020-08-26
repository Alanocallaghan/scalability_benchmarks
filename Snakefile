chains = [1, 2, 4, 8, 16, 32, 64, 128, 256]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
seeds = [7, 14, 21, 28, 35, 42]
fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# configfile: "config/snakemake_config.yaml"
# conda: config["conda"]

rule all:
    input:
        expand(
            "outputs/divide_and_conquer/d-{d}_n-{n}_s-{s}_b-{b}/",
            d = data,
            n = chains,
            b = by,
            s = seeds
        ),
        expand(
            "outputs/advi/d-{d}_s-{s}/",
            d = data,
            s = seeds
        ),
        expand(
            "outputs/downsampling/divide/d-{d}_f-{f}_s-{s}/",
            d = data,
            s = seeds,
            f = fractions
        ),
        expand(
            "outputs/downsampling/reference/d-{d}_f-{f}/",
            d = data,
            f = fractions
        ),
        expand(
            "outputs/removing/reference/d-{d}_f-{f}/",
            d = data,
            f = fractions
        ),
        expand(
            "outputs/removing/divide/d-{d}_f-{f}_s-{s}/",
            d = data,
            s = seeds,
            f = fractions
        ),
        expand(
            "outputs/batchinfo/d-{d}/",
            d = data_batch
        )



rule advi:
    output:
        "outputs/advi/d-{d}_s-{s}/"
    shell:
        """
        ./src/chain-scripts/advi.R \
            -d {wildcards.d} \
            -s {wildcards.s} \
            -o {output}
        """


rule divide_and_conquer:
    # threads: 4
    # resources:
    #     mem_mb=5000
    output:
        "outputs/divide_and_conquer/d-{d}_n-{n}_s-{s}_b-{b}/"
    shell:
        """
        ./src/chain-scripts/divide_and_conquer.R \
            -d {wildcards.d} \
            -n {wildcards.n} \
            -s {wildcards.s} \
            -b {wildcards.b} \
            -o {output}
        """

rule downsampling_ref:
    output:
        "outputs/downsampling/reference/d-{d}_f-{f}/"
    shell:
        """
        ./src/chain-scripts/downsampling_reference.R \
            -d {wildcards.d} \
            -f {wildcards.f} \
            -o {output}
        """


rule downsampling_divide:
    output:
        "outputs/downsampling/divide/d-{d}_f-{f}_s-{s}/"
    shell:
        """
        ./src/chain-scripts/downsampling_divide.R \
            -d {wildcards.d} \
            -s {wildcards.s} \
            -f {wildcards.f} \
            -o {output}
        """


rule removing_ref:
    output:
        "outputs/removing/reference/d-{d}_f-{f}/"
    shell:
        """
        ./src/chain-scripts/removing_cells_ref.R \
            -d {wildcards.d} \
            -f {wildcards.f} \
            -o {output}
        """

rule removing_divide:
    output:
        "outputs/removing/divide/d-{d}_f-{f}_s-{s}/"
    shell:
        """
        ./src/chain-scripts/removing_cells.R \
            -d {wildcards.d} \
            -s {wildcards.s} \
            -f {wildcards.f} \
            -o {output}
        """

rule batchinfo:
    output:
        "outputs/batchinfo/d-{d}/"
    shell:
        """
        ./src/chain-scripts/batchinfo.R \
            -d {wildcards.d} \
            -o {output}
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

