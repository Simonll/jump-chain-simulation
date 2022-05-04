# jump-chain simulation algorithm

Jump-chain simulation algorithm, inspired from Gillespie 1977, enables to sample subsitutional histories from complex phylogenetic models, including site-interdependent models. The timing of the substitutions (dwelltime) is sampled following a Poisson process, and the nature of the substitutions is sampled from the substitution process itself.

## How to install docker
- https://docs.docker.com/get-docker/

## Building docker containers
```bash
docker build --build-arg USER_NAME=$(whoami) --build-arg USER_ID=$(id -u ${USER}) --build-arg GROUP_ID=$(id -g ${USER}) -t ubuntu20.04/basic:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/basic --pull
docker build -t ubuntu20.04/pbmpi:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/phylobayes-mpi
docker build --build-arg CACHEBUST=$(date +%s) -t ubuntu20.04/lfp:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/LikelihoodFreePhylogenetics
docker build -t ubuntu20.04/codeml:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/codeml
docker build --build-arg USER_NAME=$(whoami) --build-arg USER_ID=$(id -u ${USER}) --build-arg GROUP_ID=$(id -g ${USER}) --build-arg CACHEBUST=$(date +%s) -t r-base3.6.3/abc:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/r-base-abc --pull
```

## How to installing snakemake
- https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

## Running snakemake
```bash
snakemake -n
snakemake -j 10 --use-conda
```
Use -s to specify snakmake file for activating different aspects of the bioinfo workflow {simu.smk, cabc.smk, codeml.smk, figures_tables.smk}. Some aspects are better run on a cluster, so some human manipulations are still expected depending of the computation environment.

## How to activate conda environment
1. from jump-chain-simulation/ create the conda environment
```bash
conda env create --file env.yml
```
2. and activate conda environment
```bash
conda activate jump_chain_simulation
```

## How to use jupyter-lab notebooks on linux based system
1. from the jump-chain-simulation/bioinfo/workflow/notebooks/ directory, launch jupyter-lab to background with:
```bash
echo "jupyter lab --no-browser --port=8888" \& | bash
```
note that port will be skip if already in use and conda environment need to be activated first \
\
2. can now launch port forwarding channel from your laptop to access the notebooks using localhost definition, noting that jupyter-lab port and localhost port forwarding need to match.
```bash
ssh -N -f -L localhost:8888:localhost:8888 user@server
```

## How to cite

BibTex Citation:
```
@article{,
    author = {Laurin-Lemay, Simon and Dickson, Kassandra and Rodrigue, Nicolas},
    title = "{Jump-chain simulation of Markov substitution processes over phylogenies}",
    journal = {},
    year = {},
    month = {},
    abstract = "{}",
    issn = {},
    doi = {},
    url = {},
    eprint = {},
}
```

## References
```
@article{Gillespie1977,
  doi = {10.1021/j100540a008},
  year = {1977},
  month = dec,
  publisher = {American Chemical Society ({ACS})},
  volume = {81},
  number = {25},
  pages = {2340--2361},
  author = {Gillespie, D T},
  title = {{Exact stochastic simulation of coupled chemical reactions}},
  journal = {The Journal of Physical Chemistry}
}
```
