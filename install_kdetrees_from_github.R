#!/usr/bin/env Rscript

# 2024-04-01
# install_kdetrees_from_github
# call to build docker image

# install /kdetrees_0.1.5.tar.gz from GitHub, as it was removed from CRAN!  
#   > Archived on 2022-05-10 as email to the maintainer is undeliverable. 
# https://cran.r-project.org/src/contrib/Archive/kdetrees/ << version 
# PV 2022-06-12s; GitHub repo 'http://github.com/grady/kdetrees'
# install.packages("kdetrees_0.1.2.tar.gz", lib='http://github.com/grady/kdetrees'  repos = NULL, type="source")


required_packages <- c("remotes")

# install the remaining required R packages
for (package in required_packages) {
  if (!require(package, character.only=T, quietly=T)) {
    sprintf("# cannot load %s, will get it from %s and install it in %s",package,repository,local_lib)
    install.packages(package, dependencies=TRUE, clean=TRUE)
  }
}


library("remotes")
#remotes::install_github("grady/kdetrees", dependencies = TRUE, force = TRUE, upgrade="always")
remotes::install_github("vinuesa/get_phylomarkers/kdetrees", dependencies = TRUE, force = TRUE, upgrade="always")

sessionInfo()
