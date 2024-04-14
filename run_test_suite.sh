#!/usr/bin/env bash

#: PROGRAM: run_test_suite.sh
#: AUTHORS: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico; @pvinmex
#: AIM: used for functional testing of the following scripts of the GET_PHYLOMARKERS package:
#: - run_get_phylomarkers_pipeline.sh
#: - estimate_pangenome_phylogenies.sh
#: - hcluster_pangenome_matrix.sh # <<< is officially distributed through the get_homologues GitHub repo

progname="${0##*/}"
version="2024-04-13"

# Activate set Bash's unofficial strict mode.
set -euo pipefail

function usage()
{
   cat <<USAGE
   
   1. Run the Makefile with:
      make
      make clean
   
   2. As standard script
   ${progname} v.${version} <full path to BASE_DIRECTORY holding core_genome and pan_genome test data>
   
   EXAMPLES:
   
     # if running from container
     ${progname} /home/you/data
   
     # if running from your host
     ${progname} $HOME/data/genomes/test_sequences

   AIM: used for functional testing of the following scripts of the GET_PHYLOMARKERS package:
        - run_get_phylomarkers_pipeline.sh
        - estimate_pangenome_phylogenies.sh
        - hcluster_pangenome_matrix.sh
   
   NOTES:
   
     Assumes that you have your test core_genome/ and pan_genome/ sequences availabe on your host machine under
      ~/data/genomes/test_sequences, or another base_dir provided as single argument to the script
  
     To let the container having acces to these data, bind mount that host directory on the container instance
      under /home/you/data with the following command to lauch the container
   
   LAUNCHING THE CONTAINER AND BIND MOUNTING THE ~/data/genomes/test_sequences dir on /home/you/data
  
     $ docker run -it --rm -v ~/data/genomes/test_sequences:/home/you/data vinuesa/get_phylomarkers:latest /bin/bash

USAGE

exit 0

}

[ $# -eq 0 ] && usage

base_dir=$1 || { echo "ERROR could not cd into $base_dir"; exit 1 ; } 

# run basic test suite for GET_PHYLOMARKERS
cd "${base_dir}"/core_genome || { echo "ERROR could not cd into $base_dir/core_genome"; exit 1 ; } 

if ls -d get_phylomarkers_run* &> /dev/null; then rm -rf get_phylomarkers_run*; fi

## NOTE: use -I 2 (but not higher, or a core may be dumped in some occasions by iqtree -T IQT_threads) instead of AUTO to speed-up tests
# 1. default on DNA sequences (uses IQ-TREE evaluating a subset of models specified in the detailed help)
echo ">>> Test #1: default run on DNA sequence ..."
run_get_phylomarkers_pipeline.sh -R 1 -t DNA -I 2 
echo

# 2. thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores and increasing k stringency 
echo ">>> Test #2 thorough FastTree searching and molecular clock analysis on DNA sequences using 10 cores and increasing k stringency ..."
run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F -k 1.2 -m 0.7 -s 20 -l 12 -T high -K -M HKY -q 0.95 
echo

# 3. test multiple models with high kdetree stringency (k=1.0) and thorogh IQT searches, using 2 seed trees
echo ">>> Test #3: multiple models with high kdetree stringency (k=1.0) and thorogh IQT searches, using 2 seed trees"
run_get_phylomarkers_pipeline.sh -R 1 -t DNA -S 'TrN,TVMe,GTR' -k 1.0 -m 0.75 -T high -N 2 -I 2
echo

# 4. IQT with proteins and moderate average bipartition support
echo ">>> Test #4: IQT with proteins and moderate average bipartition support"
run_get_phylomarkers_pipeline.sh -R 1 -t PROT -m 0.2 -I 2
echo

# 5. FastTree thorough searching on a protein dataset with thorough search
echo ">>> Test # 5. FastTree thorough searching on a protein dataset with thorough search"
run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A F -T high -m 0.2
echo

# 6. Run in population-genetics mode with K2P model, estimating the population tree with IQ-TREE
echo ">>> Test # 6. Run in population-genetics mode with K2P model, estimating the population tree on SNP matrix with IQ-TREE under best-fit model"
run_get_phylomarkers_pipeline.sh -R 2 -t DNA -S 'K2P'
echo

# 7. Run in population-genetics mode, estimating the population tree with FastTree
echo ">>> Test # 7. Run in population-genetics, estimating the population tree on SNP matrix with FastTree"
run_get_phylomarkers_pipeline.sh -R 2 -t DNA -A F
echo

# 8. estimate a ML pan-genome tree from the pan-genome matrix, using 2 independent IQT runs and UFBoot
echo ">>> Test # 8. estimate a ML pan-genome tree from the pan-genome matrix, using 2 independent IQT runs and UFBoot"
cd "${base_dir}"/pan_genome  || { echo "ERROR could not cd into $base_dir/pan_genome"; exit 1 ; }
[ -d iqtree_PGM_2_runs ] && rm -rf iqtree_PGM_2_runs
estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 2 -S UFBoot -I 2
echo

# 9. estimate a PARS  pan-genome tree with bootstrapping; 100 bootstrap replicates divided on 10 core (10 reps / core)
echo ">>> Test # 9. estimate a PARS  pan-genome tree with bootstrapping; 100 bootstrap replicates divided on 10 core (10 reps / core)"
[ -d boot_pars ] && rm -rf boot_pars
estimate_pangenome_phylogenies.sh -c PARS -R 3 -i pangenome_matrix_t0.phylip -n 10 -b 10 -j 1 -t 1
echo

# 9. cluster the pan-genome matrix and run silhoute statistic to define the optimal number of clusters
#    hcluster_pangenome_matrix.sh is officially distributed through the get_homologues GitHub repo; should not test here
#echo ">>> Test # 9. cluster the pan-genome matrix and run silhoute statistic to define the optimal number of clusters"
#hcluster_pangenome_matrix.sh -i pangenome_matrix_t0.tab -a ward.D2 -d gower -O pdf -A 'NULL,45' -X 0.8 -T "Pangenome tree"
#echo
