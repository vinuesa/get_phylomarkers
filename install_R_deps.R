#!/usr/bin/env Rscript

# check installed R packages and install only the missing ones 
# Bruno Contreras, Pablo Vinuesa, Jun2017

# Instructions to update R on Ubuntu systems, Xenial in the example:

# $ sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
# gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9
# gpg -a --export E084DAB9 | sudo apt-key add -
# $ sudo apt-get update
# $ sudo apt-get install r-base r-base-dev

# Instructions in case some of the extant packages, such as Rcpp or ape, are not up-to-date:

# $ R
# > remove.packages("Rcpp")
# > install.packages("Rcpp",dependencies=TRUE, lib="lib/R", repos="https://cloud.r-project.org")

# Instructions to install ape in MacOS in case in fails
# install gfortran from https://gcc.gnu.org/wiki/GFortranBinaries
# $ R
# > install.packages("ape",dependencies=TRUE, lib="lib/R", type="source")


repository = 'https://cloud.r-project.org'; #'http://cran.rstudio.com';
required_packages = c("ape", "kdetrees", "stringr", "vioplot", "ggplot2", "gplots", "dplyr", "seqinr")
local_lib = "./lib/R"

.libPaths( c( .libPaths(), local_lib) )

for (package in required_packages) {
  if (!require(package, character.only=T, quietly=T)) {
    sprintf("# cannot load %s, will get it from %s and install it in %s",package,repository,local_lib)
    install.packages(package, dependencies=TRUE, lib=local_lib, repos=repository)
  }
}

