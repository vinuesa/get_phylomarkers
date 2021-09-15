#!/usr/bin/env bash

progname="${0##*/}"

function usage()
{
   cat <<USAGE
   
   ${progname} <full path to BASE_DIRECTORY holding core_genome and pan_genome data computed by GET_HOMOLOGUES or similar software>
   
   EXAMPLES:
   
     # if running from container
     ${progname} /home/you/data
   
     # if running from your host
     ${progname} $HOME/data/genomes/test_sequences
   
   NOTES:
   
     Assumes that you have your test core_genome and pan_genome sequences availabe on your host machine under
      ~/data/genomes/test_sequences, or another base_dir provided as single argument to the script
  
     To let the container having acces to these data, bind mount that host directory on the container instance
      under /home/you/data with the following command to lauch the container
   
   LAUNCHING THE CONTAINER AND BIND MOUNTING THE ~/data/genomes/test_sequences dir on /home/you/data
  
     $ docker run -it --rm -v ~/data/genomes/test_sequences:/home/you/data get_phylomarkers:latest /bin/bash

USAGE

exit 0

}

[ $# -eq 0 ] && usage

base_dir=$1 || { echo "ERROR could not cd into $base_dir"; exit 1 ; } 

# run basic test suite for GET_PHYLOMARKERS
cd "${base_dir}"/core_genome || { echo "ERROR could not cd into $base_dir/core_genome"; exit 1 ; } 
run_get_phylomarkers_pipeline.sh -R 1 -t DNA -k 1.0
run_get_phylomarkers_pipeline.sh -R 1 -t DNA -A F
run_get_phylomarkers_pipeline.sh -R 2 -t DNA

run_get_phylomarkers_pipeline.sh -R 1 -t PROT
run_get_phylomarkers_pipeline.sh -R 1 -t PROT -A F

cd "${base_dir}"/pan_genome  || { echo "ERROR could not cd into $base_dir/pan_genome"; exit 1 ; }
estimate_pangenome_phylogenies.sh -f pangenome_matrix_t0.fasta -r 1 -S UFBoot
