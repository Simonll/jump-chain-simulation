# jump-chain-simulation

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
------------------------
