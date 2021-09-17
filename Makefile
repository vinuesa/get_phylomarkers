# Makefile for Travis CI tests and installation
# https://docs.travis-ci.com/user/languages/perl

test:
	perl test_get_phylomarkers.t

install:
	Rscript install_R_deps.R

clean:
	rm -rf test_sequences/core_genome/get_phylomarkers_run* 
	rm -rf test_sequences/pan_genome/iqtree_PGM_*
	rm -rf test_sequences/pan_genome/boot_pars
