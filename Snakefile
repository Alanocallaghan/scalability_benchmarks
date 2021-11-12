chains = [1, 2, 4, 8, 16, 32, 64]
chains_timing = [2, 4, 8, 16, 32, 64, 128]
by = ["gene"]
data = ["buettner", "chen", "tung", "zeisel"]
data_all = ["buettner", "chen", "tung", "zeisel", "ibarra-soria"]
data_spikes = ["buettner", "tung", "zeisel"]
data_batch = ["tung", "zeisel"]
seeds = [14, 21, 28, 35, 42]
# seeds = [42]
fractions = [x/10 for x in range(2, 11, 2)]
iterations = 40000


shell.prefix("source src/modules.sh; ")

rule all:
    input:
        "figs/diffexp_plot.pdf",
        "figs/overlap_diff_genes.pdf",
        "figs/libsize_density.pdf",
        "figs/complexity_density.pdf",
        "figs/expression_density.pdf",
        "figs/dropout_density.pdf",
        "figs/true_positives.pdf",
        "figs/time_plot.pdf",
        # "figs/de_batch_tung.pdf",
        # "figs/de_batch_zeisel.pdf",
        # "figs/ess_batch.pdf",
        "figs/removing_cells.pdf",
        "figs/downsampling.pdf",
        "figs/elbo/tung.pdf",
        "figs/elbo/buettner.pdf",
        "figs/elbo/zeisel.pdf",
        "figs/elbo/chen.pdf",
        "figs/hpd_mu.pdf",
        "figs/fixnu-comparison.pdf",
        "figs/hpd_epsilon.pdf",
        "tables/data-summary.tex",
        expand(
            "figs/hpd/{dataset}",
            dataset = data
        )
        # ,
        # "figs/ess_mu.pdf",
        # "figs/ess_epsilon.pdf",
        
        # "figs/de_id_tung.pdf",
        # "figs/de_id_zeisel.pdf",
        # "figs/ess_id.pdf",


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
            dataset = data,
            seed  = seeds
        ),
        expand(
            "outputs/downsampling/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/",
            dataset = "zeisel",
            seed = seeds,
            fraction = fractions
        ),
        expand(
            "outputs/downsampling/reference/data-{dataset}_fraction-{fraction}/",
            dataset = "zeisel",
            fraction = fractions
        ),
        expand(
            "outputs/removing/reference/data-{dataset}_fraction-{fraction}/",
            dataset = "zeisel",
            fraction = fractions
        ),
        expand(
            "outputs/removing/divide/data-{dataset}_fraction-{fraction}_seed-{seed}/",
            dataset = "zeisel",
            seed = seeds,
            fraction = fractions
        ),
        # expand(
        #     "outputs/batchinfo/data-{dataset}/",
        #     dataset = data_batch
        # ),
        expand(
            "outputs/true-positives/data-ibarra-soria_nsubsets-{nsubsets}_seed-{seed}.rds",
            nsubsets = chains,
            seed = seeds
        ),
        expand(
            "outputs/time/{dataset}_{n}.rds",
            dataset = data_all,
            n = chains_timing
        )
    output:
        "figs/diffexp_plot.pdf",
        "figs/overlap_diff_genes.pdf",
        "figs/norm_plot.pdf",
        "figs/hpd_mu.pdf",
        "figs/hpd_epsilon.pdf",
        "rdata/main_done.RData"
    shell:
        """
        Rscript src/analysis/main.R
        """

# rule identifiability_plot:
#     input: 
#     output:
#         "figs/de_id_tung.pdf",
#         "figs/de_id_zeisel.pdf",
#         "figs/ess_id.pdf"
#     shell:
#         """
#         Rscript src/analysis/identifiability.R
#         """

rule batchinfo_plot:
    input:
        expand(
            "outputs/batchinfo/data-{dataset}/",
            dataset = data_batch
        )
    output:
        "figs/de_batch_tung.pdf",
        "figs/de_batch_zeisel.pdf",
        "figs/ess_batch.pdf"
    shell:
        """
        Rscript src/analysis/batchinfo.R
        """

rule true_positive_plot:
    input:
        dc = expand(
            "outputs/true-positives/data-ibarra-soria_nsubsets-{nsubsets}_seed-{seed}.rds",
            nsubsets = chains,
            seed = seeds
        ),
        advi = expand(
            "outputs/true-positives/data-ibarra-soria_advi-{seed}.rds",
            seed = seeds
        )
    output: 
        "figs/true_positives.pdf"
    shell:
        """
        Rscript src/analysis/true_positives.R
        """

# rule hpd_ess_plot:
#     input:
#         expand(
#             "outputs/divide_and_conquer/data-{dataset}_nsubsets-{nsubsets}_seed-{seed}_by-{by}/",
#             dataset = data,
#             nsubsets = chains,
#             by = by,
#             seed = seeds
#         )
#     output: 
#         expand(
#             "figs/{metric}_{parameter}.pdf",
#             metric = ["hpd", "ess"],
#             parameter = ["mu", "epsilon"]
#         )
#     shell:
#         """
#         Rscript src/analysis/hpd.R
#         Rscript src/analysis/ess.R
#         """

rule hpd_comparison_plot:
    resources: mem_mb=20000
    input:
        reference = "outputs/divide_and_conquer/data-{dataset}_nsubsets-16_seed-14_by-gene/chains.rds",
        divide = "outputs/divide_and_conquer/data-{dataset}_nsubsets-1_seed-14_by-gene/chains.rds",
        advi = "outputs/advi/data-{dataset}_seed-14/chain.rds"
    output:
        directory("figs/hpd/{dataset}")
    shell:
        """
        Rscript src/analysis/hpd_comparison.R \
            -d {wildcards.dataset} \
            -o {output}
        """

rule removing_cells_plot:
    resources: mem_mb=20000
    input:
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
        )
    output:
        "figs/removing_cells.pdf"
    shell:
        """
        Rscript src/analysis/removing_cells.R
        """

rule downsampling_plot:
    input:
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
        )
    output: 
        "figs/downsampling.pdf"
    shell: 
        """
        Rscript src/analysis/downsampling.R
        """


rule time_plot:
    input:
        expand(
            "outputs/advi/data-{dataset}_seed-{seed}/",
            dataset = data_spikes,
            seed  = seeds
        ),
        expand(
            "outputs/time/{dataset}_{n}.rds",
            dataset = data,
            n = chains_timing
        )
    output:
        "figs/time_plot.pdf"
    shell:
        """
        Rscript src/analysis/time_plot.R
        """

rule data_comparison:
    input:
        expand(
            "rdata/{data}.rds",
            data = data
        )
    output:
        "figs/dropout_density.pdf",
        "figs/expression_density.pdf",
        "figs/complexity_density.pdf",
        "figs/libsize_density.pdf"
    shell:
        """
        Rscript src/analysis/data_comparison.R
        """


rule data_summary:
    resources: mem_mb=20000
    output:
        "tables/data-summary.tex"
    input:
        "rdata/zeisel.rds",
        "rdata/chen.rds",
        "rdata/buettner.rds",
        "rdata/ibarra-soria.rds",
        "rdata/tung.rds"
    shell:
        """
        Rscript ./src/analysis/data_summary.R
        """

rule advi:
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
    resources: mem_mb=15000, runtime=10000
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


rule true_positives_advi:
    resources: mem_mb=20000, runtime=5000
    input:
        "rdata/ibarra-soria.rds"
    output:
        "outputs/true-positives/data-ibarra-soria_advi-{seed}.rds"
    shell:
        """
        Rscript ./src/chains/true_positives_advi.R \
            --seed {wildcards.seed} \
            --output {output}
        """


rule batchinfo:
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


# rule cell:
#     resources: mem_mb=10000, runtime=3000
#     input:
#         "rdata/{dataset}.rds"
#     output:
#         "outputs/cell_splitting/{dataset}.rds"
#     shell:
#         """
#         Rscript ./src/chains/batchinfo.R \
# 	        --iterations {iterations} \
#             --data {wildcards.dataset} \
#             --output {output}
#         """


rule plot_fixnu:
    resources: mem_mb=10000, runtime=3000
    input:
        "outputs/fix_nu/"
    output:
        "figs/fixnu-comparison.pdf"
    shell:
        """
        Rscript ./src/analysis/plot_fix_nu.R -i {input}
        """


rule fixnu:
    resources: mem_mb=10000, runtime=3000
    input:
        "rdata/chen.rds"
    output:
        directory("outputs/fix_nu/")
    shell:
        """
        Rscript ./src/chains/fix_nu.R \
	        --iterations {iterations} \
            --output {output}
        """


rule data:
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

rule elbo_plots:
    input:
        expand(
            "outputs/advi/data-{dataset}_seed-{seed}/",
            seed = seeds,
            allow_missing = True
        )
    output:
        "figs/elbo/{dataset}.pdf"
    shell:
        """
        Rscript src/analysis/elbo_plots.R
        """
