#!/usr/bin/env Rscript

# check R packages required by get_phylomarkers and install missing ones 
# Bruno Contreras Moreira, Pablo Vinuesa, Nov2017-Oct2020
# version: 2024-04-15

# WARNING: some packages require C (gcc) and C++ (g++) compilers to be installed
# These can be installed with these commands:
# sudo apt-get install g++      # Ubuntu
# sudo yum install gcc-c++      # CentOS   

# Instructions to update R on Ubuntu systems, bionic-cran35 in the example:

# $ sudo echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" | sudo tee -a /etc/apt/sources.list
# gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
# gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -
# $ sudo apt-get update
# $ sudo apt-get install r-base r-base-dev

# Instructions in case some preinstalled packages, such as Rcpp or ape, are not up-to-date:

# $ R
# > remove.packages("Rcpp")
# > install.packages("Rcpp",dependencies=TRUE, lib="lib/R", repos="https://cloud.r-project.org")

# Instructions to install ape from source, in case it causes kdetrees errors
# In MacOS this requires installing gfortran from https://gcc.gnu.org/wiki/GFortranBinaries
# $ R
# > install.packages("ape", dependencies=TRUE, lib="lib/R", type="source")
# > install.packages("kdetrees",dependencies=TRUE, lib="lib/R", type="source")

repository = 'https://CRAN.R-project.org' # ; 'http://cran.rstudio.com' # 'https://cloud.r-project.org'; https://CRAN.R-project.org/

# do not change (reduce), as it includes dependencies for the GET_HOM+GET_PHYLO image,
# and from v2.3.0 (2021-09-18) the GET_PHYLO package also includes hcluster_pangenome_matrix.sh, which require "cluster", "dendextend", "factoextra"
# Note that plyr should be called before dplyr; stringi before stringr;
# https://github.com/tidyverse/stringr/issues/320
# 
# required_packages = c("devtools", "ape", "cluster", "gplots", "vioplot", "plyr", "dplyr", "ggplot2", "stringi", "stringr", "seqinr", "dendextend", "factoextra")
# 2024-04-13: removed "cluster", "dendextend", "factoextra", which are only required by hcluster_pangenome_matrix.sh, 
#    which is distributed through the GET_HOMOLOGUES GitHub repo.
# 2024-04-15: replaced devtools for the much lighter remotes.
required_packages = c("remotes", "ape", "gplots", "vioplot", "plyr", "dplyr", "ggplot2", "stringi" "stringr", "seqinr")


local_lib = "./lib/R"

# Note: this command should be added to .Rprofile to make permanent, by calling it with each new shell/session start
.libPaths( c( .libPaths(), local_lib) )

# make sure we get latest ape Rcpp packages installed
#remove.packages(c("ape", "cluster", "dendextend", "factoextra", "kdetrees"), lib=local_lib)

# Install ape && kdetrees from source
# install.packages(c("ape", "kdetrees"), dependencies=TRUE, lib="lib/R", type="source")

# install the remaining required R packages
for (package in required_packages) {
  if (!require(package, character.only=T, quietly=T)) {
    sprintf("# cannot load %s, will get it from %s and install it in %s",package,repository,local_lib)
    install.packages(package, dependencies=TRUE, lib=local_lib, repos=repository, clean=TRUE)
  }
}

# install /kdetrees_0.1.5.tar.gz from GitHub, as it was removed from CRAN!  
#   > Archived on 2022-05-10 as email to the maintainer is undeliverable. 
# https://cran.r-project.org/src/contrib/Archive/kdetrees/ << version 
# PV 2022-06-12s; GitHub repo 'http://github.com/grady/kdetrees'
# install.packages("kdetrees_0.1.5.tar.gz", lib='/usr/lib/R/site-library', repos = NULL, type="source")

library("devtools")
install_github("grady/kdetrees", dependencies = TRUE)

sessionInfo()
