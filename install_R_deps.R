#!/usr/bin/Rscript

# checks installed R packages and installs only the missing ones 

repository = 'http://cran.rstudio.com/'; #'https://cloud.r-project.org'

required_packages = c("ape", "kdetrees", "stingr", "vioplot", "ggplot2", "gplots", "plyr", "seqinr")
new.packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE, lib="lib/R", repos=repository)
