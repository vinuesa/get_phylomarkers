# version 2021-09-09
# use ubuntu:18.04 instead of ubuntu:20.04 to avoid problems with R's kdetree package
FROM ubuntu:18.04
FROM r-base:3.6.3

LABEL maintainer="vinuesa[at]ccg dot unam dot mx"
LABEL version="0.2_2021-09-09"
LABEL description="Image for GET_PHYLOMARKERS: a software package designed to identify optimal genomic markers for phylogenomics, \
population genetics and genomic taxonomy. It implements a pipeline to filter orthologous gene clusters computed by the companion package \
GET_HOMOLOGUES to select those with optimal phylogenetic attributes. Top-scoring alignments are concatenated into a supermatrix, \
which is used to estimate the species tree under the maximum-likelihood (ML) criterion with state-of-the-art fast ML tree searching algorithms. \
GET_PHYLOMARKERS can also estimate ML and parsimony trees from the pan-genome matrix, including unsupervised learning methods to determine the \
optimal number of clusters from pan-genome and average genomic distance matrices. A detailed manual and step-by-step tutorials document the software \
and help the user to get quickly up and running. For your convenience, html and markdown versions of the documentation material are available."


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

# clone get_phylomarkers and install required R packages (slow)
# Don't forget to update to latest github cachebuster=b953b30
RUN cachebuster=3853f87 git clone https://github.com/vinuesa/get_phylomarkers.git
RUN cd get_phylomarkers && Rscript install_R_deps.R
RUN cd get_phylomarkers && git pull

# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare user env
RUN useradd --create-home --shell /bin/bash you && usermod -aG docker you
# set USER=you for pgrep calls, as the only ENV-VARS documented to be set are HOME HOSTNAME PATH and TERM
#  https://docs.docker.com/engine/reference/run/#env-environment-variables
ENV USER=you
USER you

WORKDIR /home/you
ENV PATH="/get_phylomarkers:${PATH}"


