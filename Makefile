# Makefile for Travis CI tests and installation
# https://docs.travis-ci.com/user/languages/perl

test:
	perl test_get_phylomarkers.t
	rm -rf test_sequences/core_genome/get_phylomarkers_run* 

install:
	Rscript install_R_deps.R
