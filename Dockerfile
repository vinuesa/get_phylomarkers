# Dockerfile version 2021-09-13
# - build images using as context the freshly pulled git repo directory get_phylomarkers

## Base images
# Use ubuntu:18.04 as base layer and r-base:3.6.3-bionic to avoid problems with R's kdetree package
FROM ubuntu:18.04

## Add metadata 
LABEL   maintainer="Pablo Vinuesa @pvinmex <https://www.ccg.unam.mx/~vinuesa/>" \
    software="GET_PHYLOMARKERS" \
    version="v2.2.9.3_13sep21" \
    note="This version used .dockerignore to exclude several src dirs, including the test_sequences. \
    Build environment is freshly cloned GitHub repo <https://github.com/vinuesa/get_phylomarkers>." \
    about.summary="an open source tool to estimate maximum-likelihood core-genome phylogenies and pan-genome trees" \
    about.home="https://github.com/vinuesa/get_phylomarkers" \
    about.documentation="https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual" \
    about.license="GPL-3.0" \
    about.license_file="https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE" \
    about.tags="Phylogenetics"

## Install required linux tools
RUN apt-get update && apt-get install --no-install-recommends -y \
apt-utils \
bash-completion \
bc \
cpanminus \
curl \
htop \
libprocps6 \
procps \
wget \
&& apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && cpanm Term::ReadLine

## use r-base:3.6.3-bionic to avoid problems with R's kdetree package
FROM rstudio/r-base:3.6.3-bionic

## mkdir get_phylomarkes in /, copy all contents into it; make it the working directory and install required R packages
RUN mkdir get_phylomarkers 
COPY . /get_phylomarkers 
WORKDIR /get_phylomarkers 
RUN Rscript /get_phylomarkers/install_R_deps.R

## add version tag to image
ARG version
LABEL version=$version
RUN echo $version

## Prepare USER env
RUN useradd --create-home --shell /bin/bash you && usermod -aG root you
# set USER=you for pgrep calls, as the only ENV-VARS documented to be set are HOME HOSTNAME PATH and TERM
#  https://docs.docker.com/engine/reference/run/#env-environment-variables
ENV USER=you
USER you

WORKDIR /home/you
ENV PATH="${PATH}:/get_phylomarkers"

# make sure the user gets a Bash shell and some help to get started
CMD ["/bin/bash"]

## USAGE:
# 1. Copy test data from /get_phylomarkers/test_sequences
# - open a terminal on your post and type: 
# $ cd && mkdir -p ~/data/genomes/test_sequences

# 2. Assuming that you have your core_genome and/or pan_genome data availabe in ~/data/genomes/test_sequences 
# bind mount that host directory on the container instance under /home/you/data with the following command:
# $ docker run -it --rm -v ~/data/genomes/test_sequences:/home/you/data get_phylomarkers:latest /bin/bash

# Copy the test sequences to your data directory and extract them from the gzipped tar files
# $ cp /get_phylomarkers/test_sequences/*tgz ~/data/
# $ cd ~/data 
# $ for f in *tgz; do tar -xzf $f; done
# $ rm *tgz && ls

## Functional testing
# $ cd ~/data/core_genome
# $ ls
# $ run_get_phylomarkers_pipeline.sh -h
# $ run_get_phylomarkers_pipeline.sh -R 1 -t DNA
# $ cd ~/data/pan_genome
# $ estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 1 -S UFBoot

## Functional testing can also be run as follows; 
# Note: if you ran the previous block, you will need to remove the get_phylomarkers* dir in core_genome 
#       and iqtree* dir in pan_genome
# $ run_test_suite.sh /home/you/data

