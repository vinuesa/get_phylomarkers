# version 2021-09-09
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
  wget \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get -y autoremove

# clone get_phylomarkers and install required R packages (slow)
#RUN cachebuster=b953b30 git clone https://github.com/vinuesa/get_phylomarkers.git
RUN cachebuster=ea8896d git clone https://github.com/vinuesa/get_phylomarkers.git
RUN cd get_phylomarkers && Rscript install_R_deps.R
RUN cd get_phylomarkers && git pull

# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare user env
RUN useradd -ms /bin/bash you
USER you
WORKDIR /home/you
ENV PATH="/get_phylomarkers:${PATH}"


