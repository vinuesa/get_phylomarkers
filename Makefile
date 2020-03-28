# Makefile for Travis CI tests and installation
# https://docs.travis-ci.com/user/languages/perl

test:
	perl test_get_phylomarkers.t
	#rm -rf sample_plasmids_gbk_homologues/tmp/ 
