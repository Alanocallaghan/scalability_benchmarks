# Benchmarking of scalable Bayesian inference methods for BASiCS

This is the repo used for all of the benchmarking of ADVI, divide and conquer,
and, in a limited way, HMC.
It's located at:
`/exports/igmm/eddie/cvallejo-scRNAseq/alan/scalability_benchmarks/`

# Setup

First install miniconda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Then make a conda environment containing snakemake for later.
```
conda create -n snakemake snakemake
```

If the following conda environment is available, use it:
```
/exports/igmm/eddie/cvallejo-scRNAseq/alan/conda/envs/scalability
```

If not, create the conda environment used in the analyses:

```
source ./src/mkconda.sh
```

# Output directory

Most of the outputs are written to the creatively-named `outputs` directory. I
have this set up to be in my scratch space, because the ouputs are large.

```
mkdir -p /exports/eddie/scratch/`whoami`/scalability
ln -s /exports/eddie/scratch/`whoami`/scalability outputs
```

Graphical and tabular outputs are saved to `tables` and `figs`.

# Running

Running scripts is handled by snakemake. I use 
[a profile](https://github.com/Alanocallaghan/snakemake-sync-bq-sub)
to handle resource allocation etc.

Then I'd just `ssh wildwest`.
To actually run jobs, you need to specify the job limit and the profile:

```
conda activate snakemake
snakemake --jobs 250 --profile cluster-sync
```


# eddie-specific

On eddie you first you need to connect to a wild west node. I use `node2c16`.
I have this set up in my `~/.ssh/config` as follows:

```
Host eddie
    HostName eddie.ecdf.ed.ac.uk
    User s1372510
    ForwardX11 yes

Host wildwest
    Hostname node2c16
    User s1372510
    ProxyJump eddie
    ForwardX11 yes
```
