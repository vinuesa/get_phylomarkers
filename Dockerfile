# version 2021-09-11
# build using as context the freshly pulled git repo directory get_phylomarkers
#   * this seems to be the only way to make .dockerignore actually ingore the .git dir
# use ubuntu:18.04 instead of ubuntu:20.04 to avoid problems with R's kdetree package
## Start multistage build
# Stage 0
FROM ubuntu:18.04 AS ubuntu_18.04

LABEL maintainer="Pablo Vinuesa <vinuesa[at]ccg dot unam dot mx>"
LABEL software="GET_PHYLOMARKERS"
LABEL version="0.3_2021-09-11"
LABEL note="This version used .dockerignore to exclude several src dirs, including the test_sequences. \
 Build environment is freshly cloned GitHub repo."
LABEL about.summary="an open source tool to estimate maximum-likelihood core-genome phylogenies and pan-genome trees"
LABEL about.description="Image for GET_PHYLOMARKERS: a software package designed to identify optimal genomic markers for phylogenomics, \
population genetics and genomic taxonomy. It implements a pipeline to filter orthologous gene clusters computed by the companion package \
GET_HOMOLOGUES to select those with optimal phylogenetic attributes. Top-scoring alignments are concatenated into a supermatrix, \
which is used to estimate the species tree under the maximum-likelihood (ML) criterion with state-of-the-art fast ML tree searching algorithms. \
GET_PHYLOMARKERS can also estimate ML and parsimony trees from the pan-genome matrix, including unsupervised learning methods to determine the \
optimal number of clusters from pan-genome and average genomic distance matrices. A detailed manual and step-by-step tutorials document the software \
and help the user to get quickly up and running. For your convenience, html and markdown versions of the documentation material are available."
LABEL about.home="https://github.com/vinuesa/get_phylomarkers"
LABEL about.documentation=""
LABEL about.license="GPL-3.0"
LABEL about.license_file="https://github.com/vinuesa/get_phylomarkers/blob/master/LICENSE"
LABEL about.tags="Bioinformatics Phylogenetics Genomics Pan-genomics Ubuntu"


# Install dependencies from repos: GD for graphics, libidn11 for BLAST+
RUN apt-get update && apt-get install -y --no-install-recommends \
  apt-utils \
  bc \
  curl \
  git \
  htop \
  procps \
  wget \
  && rm -rf /var/lib/apt/lists/*

## Multistage build, stage 1
FROM r-base:3.6.3 AS r-base3.6.3

## mkdir get_phylomarkes in /, copy all contents into it and make it the working directory
RUN mkdir get_phylomarkers 
COPY . /get_phylomarkers
WORKDIR /get_phylomarkers
RUN Rscript /get_phylomarkers/install_R_deps.R

## add version name to image
ARG version
LABEL version=$version
RUN echo $version

## Prepare USER env
RUN useradd --create-home --shell /bin/bash you && usermod -aG docker you
# set USER=you for pgrep calls, as the only ENV-VARS documented to be set are HOME HOSTNAME PATH and TERM
#  https://docs.docker.com/engine/reference/run/#env-environment-variables
ENV USER=you
USER you

WORKDIR /home/you
ENV PATH="/get_phylomarkers:${PATH}"

# make sure the user gets a Bash shell
CMD ["/bin/bash"]

