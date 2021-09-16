## Dockerfile version 2021-09-15 <vivaMX>
# - build images using as context the freshly pulled git repo directory get_phylomarkers

## Base images
# Use ubuntu:18.04 as base layer and r-base:3.6.3-bionic to avoid problems with R's kdetree package
#  Note however, that r-base:3.6.3-bionic overwrites python2.7, which needs to be installed after the former
FROM ubuntu:18.04

## Add metadata 
LABEL maintainer="Pablo Vinuesa <vinuesa@ccg.unam.mx>" \
 software="get_phylomarkers" \
 version="20210915" \
 summary="An open source tool to estimate maximum-likelihood core-genome phylogenies and pan-genome trees" \
 home="https://github.com/vinuesa/get_phylomarkers" \
 about.documentation="https://vinuesa.github.io/get_phylomarkers/#get_phylomarkers-manual" \
 about.license="GPL-3.0" \
 about.license_file="https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE" \
 about.tags="bioinformatics genomics phylogenetics pipeline ubuntu"

## Install required linux tools
RUN apt-get update && apt-get install --no-install-recommends -y \
bash-completion \
bc \
build-essential \
cpanminus \
gcc \
make \
wget \
&& apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && cpanm Term::ReadLine

## use r-base:3.6.3-bionic to avoid problems with R's kdetree package
# NOTE: seems that rstudio/r-base:3.6.3-bionic overwritesÂpython2.7
FROM rstudio/r-base:3.6.3-bionic

## mkdir get_phylomarkes in /, copy all contents into it; make it the working directory & install required R packages
RUN mkdir get_phylomarkers 
COPY . /get_phylomarkers 
WORKDIR /get_phylomarkers 
RUN Rscript /get_phylomarkers/install_R_deps.R

## python2.7 required by paup; python2.7-dev to get libpython2.7.so.1.0
#   add python2.7 at this stage, as rstudio/r-base:3.6.3-bionic seems to overwriteÂpython2.7 in the first RUN apt-get above
#   this seems the only way to get python2.7 and python2.7-dev in newer ubuntu:18.04 images
#   as otherwise only python 3.6 gets installed.
#  Copy libnw required by newick-utils to /usr/local/lib and export LD_LIBRARY_PATH for ldconfig
RUN apt-get update && apt-get install --no-install-recommends -y python2.7 python2.7-dev \
&& apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
&& cp /get_phylomarkers/lib/libnw.so /usr/local/lib \
&& export LD_LIBRARY_PATH=/usr/local/lib \
&& ldconfig

## add version tag to image
ARG version
LABEL version=$version
RUN echo $version

## Prepare USER env
RUN useradd --create-home --shell /bin/bash you
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
# $ docker run -it --rm -v ~/data/genomes/test_sequences:/home/you/data vinuesa/get_phylomarkers:latest /bin/bash

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

## Thorough functional testing (8 tests calling the two main scripts) can also be run automatically as follows,
#   assuming that you have installed the core_genome and pan_genome test sequences, as indicated above.
# $ cd && run_test_suite.sh /home/you/data

