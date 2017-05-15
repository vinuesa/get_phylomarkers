#!/usr/bin/Rscript

# checks installed R packages and installs only the missing ones 

# Instructions to update R on Ubuntu systems, Xenial in the example:

# $ sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
# gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9
# gpg -a --export E084DAB9 | sudo apt-key add -
# $ sudo apt-get update
# $ sudo apt-get install r-base r-base-dev

# Instructions in case some of the extant packages are not up-to-date:

# $ R
# > remove.packages("Rcpp")
# > install.packages("Rcpp",dependencies=TRUE, lib="lib/R", repos="https://cloud.r-project.org")

repository = 'https://cloud.r-project.org'; #'http://cran.rstudio.com/';

required_packages = c("ape", "kdetrees", "stringr", "vioplot", "ggplot2", "gplots", "plyr", "seqinr")
new.packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE, lib="lib/R", repos=repository)
