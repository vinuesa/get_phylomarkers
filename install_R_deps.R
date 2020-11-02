#!/usr/bin/env Rscript

# check R packages required by get_phylomarkers and install missing ones 
# Bruno Contreras Moreira, Pablo Vinuesa, Nov2017-Oct2020

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
# > install.packages("ape",dependencies=TRUE, lib="lib/R", type="source")
# > install.packages("kdetrees",dependencies=TRUE, lib="lib/R", type="source")

repository = 'https://cloud.r-project.org'; #'http://cran.rstudio.com'

required_packages = c("ape", "kdetrees", "stringr", "vioplot", "ggplot2", "gplots", "dplyr", "seqinr", "plyr")

local_lib = "./lib/R"

.libPaths( c( .libPaths(), local_lib) )

for (package in required_packages) {
  if (!require(package, character.only=T, quietly=T)) {
    sprintf("# cannot load %s, will get it from %s and install it in %s",package,repository,local_lib)
    install.packages(package, dependencies=TRUE, lib=local_lib, repos=repository)
  }
}

sessionInfo()
